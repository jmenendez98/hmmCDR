import argparse
import ast
import concurrent.futures
import os

import numpy as np
from hmmlearn import hmm

from hmmCDR.bed_parser import bed_parser
from hmmCDR.calculate_matrices import calculate_matrices


class hmmCDR:
    def __init__(
        self,
        n_iter,
        tol,
        min_cdr_size,
        min_cdr_score,
        min_low_conf_score,
        main_color,
        low_conf_color,
        output_label,
    ):

        self.n_iter = n_iter
        self.tol = tol

        self.min_cdr_size = min_cdr_size
        self.min_cdr_score = min_cdr_score
        self.min_low_conf_score = min_low_conf_score

        self.main_color = main_color
        self.low_conf_color = low_conf_color
        self.output_label = output_label

    def runHMM(self, methylation, e_mtx, t_mtx):
        # create HMM
        model = hmm.CategoricalHMM(n_components=2, n_iter=self.n_iter, tol=self.tol, init_params="")

        # initialize HMM variables
        model.startprob_ = np.array([1.0, 0.0])
        model.emissionprob_ = e_mtx
        model.transmat_ = t_mtx

        # get HMM scoring from running HMM on data
        emission_data = methylation["emissions"].reshape(-1, 1)
        _, predicted_states = model.decode(emission_data, algorithm="viterbi")
        log_likelihood, responsibilities = model.score_samples(emission_data)

        methylation["hmm_cdr_score"] = np.array([score[1]*100 for score in responsibilities], dtype=float) 

        # no call is 0, low confidence is 1, high confidence is 2
        methylation["cdr_calls"] = np.zeros(len(methylation["hmm_cdr_score"]), dtype=int) 
        low_conf_idx = np.where(methylation["hmm_cdr_score"] >= self.min_low_conf_score)[0]
        methylation["cdr_calls"][low_conf_idx] = 1
        high_conf_idx = np.where(methylation["hmm_cdr_score"] >= self.min_cdr_score)[0]
        methylation["cdr_calls"][high_conf_idx] = 2

        return methylation

    def call_hmm_cdrs(self, methylation_emissions_priors_hmm):
        cdr_entries = {"starts": [], "ends": [], "name": [], "score": [], "strand": [], "itemRgb": []}

        methyl_pos = np.array(methylation_emissions_priors_hmm["starts"], dtype=int)
        methyl_hmm_scores = np.array(methylation_emissions_priors_hmm["hmm_cdr_score"], dtype=float)
        cdr_calls = np.array(methylation_emissions_priors_hmm["cdr_calls"], dtype=int)

        high_conf_cdrs_idx = np.where(cdr_calls == 2)[0]
        high_conf_cdrs_idx_diff = np.diff(high_conf_cdrs_idx)
        high_conf_cdrs_idx_breaks = np.where(high_conf_cdrs_idx_diff != 1)[0] + 1 
        high_conf_runs = np.split(high_conf_cdrs_idx, high_conf_cdrs_idx_breaks)

        for run in high_conf_runs:
            if len(run) == 0:
                continue

            min_idx, max_idx = np.min(run), np.max(run)
            start, end = methyl_pos[min_idx], methyl_pos[max_idx]
            run_length = end - start

            slice_scores = methyl_hmm_scores[min_idx:max_idx+1]
            if slice_scores.size == 0:
                continue
            score = np.mean(slice_scores)
            if np.isnan(score):
                continue

            if run_length >= self.min_cdr_size:
                cdr_entries["starts"].append(start)
                cdr_entries["ends"].append(end)
                cdr_entries["name"].append(f"sub{self.output_label}")
                cdr_entries["strand"].append(".")
                cdr_entries["score"].append(score)
                cdr_entries["itemRgb"].append(self.main_color)
            else:
                cdr_entries["starts"].append(methyl_pos[min_idx-1]+1)
                cdr_entries["ends"].append(methyl_pos[min_idx+1]-1)
                cdr_entries["name"].append(f"low_conf_sub{self.output_label}")
                cdr_entries["strand"].append(".")
                cdr_entries["score"].append(score)
                cdr_entries["itemRgb"].append(self.low_conf_color)

            low_conf_cdrs_idx = np.where(cdr_calls == 1)[0]
            low_conf_cdrs_idx_diff = np.diff(low_conf_cdrs_idx)
            low_conf_cdrs_idx_breaks = np.where(low_conf_cdrs_idx_diff != 1)[0] + 1 
            low_conf_runs = np.split(low_conf_cdrs_idx, low_conf_cdrs_idx_breaks)

        for run in low_conf_runs:
            if len(run) == 0:
                continue

            min_idx, max_idx = np.min(run), np.max(run)
            start, end = methyl_pos[min_idx], methyl_pos[max_idx]
            run_length = end - start

            slice_scores = methyl_hmm_scores[min_idx:max_idx+1]
            if slice_scores.size == 0:
                continue
            score = np.mean(slice_scores)
            if np.isnan(score):
                continue

            cdr_entries["starts"].append(methyl_pos[min_idx-1]+1)
            cdr_entries["ends"].append(methyl_pos[min_idx+1]-1)
            cdr_entries["name"].append(f"low_conf_sub{self.output_label}")
            cdr_entries["strand"].append(".")
            cdr_entries["score"].append(score)
            cdr_entries["itemRgb"].append(self.low_conf_color)

        return cdr_entries

    def hmm_single_chromosome(self, chrom, methylation_emissions_priors, emission_matrix, transition_matrix):
        methylation_emissions_priors_hmm = self.runHMM(methylation_emissions_priors, 
                                                       e_mtx=emission_matrix, 
                                                       t_mtx=transition_matrix)
        cdrs = self.call_hmm_cdrs(methylation_emissions_priors_hmm)

        return chrom, cdrs, methylation_emissions_priors_hmm

    def hmm_all_chromosomes(self, methylation_emissions_priors_all_chroms, emission_matrix_all_chroms, transition_matrix_all_chroms):
        cdrs_all_chroms = {}
        methylation_emissions_priors_hmm_all_chroms = {}

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.hmm_single_chromosome,
                    chrom,
                    methylation_emissions_priors_all_chroms[chrom],
                    emission_matrix_all_chroms[chrom],
                    transition_matrix_all_chroms[chrom],
                ): chrom
                for chrom in methylation_emissions_priors_all_chroms
            }

            for future in concurrent.futures.as_completed(futures):
                (
                    chrom,
                    cdrs,
                    methylation_emissions_priors_hmm,
                ) = future.result()

                cdrs_all_chroms[chrom] = cdrs
                methylation_emissions_priors_hmm_all_chroms[chrom] = methylation_emissions_priors_hmm

        return cdrs_all_chroms, methylation_emissions_priors_hmm_all_chroms


def parse_command_line_arguments():
    argparser = argparse.ArgumentParser(
        description="Process input files with optional parameters."
    )

    argparser.add_argument("bedmethyl", type=str, help="Path to the bedMethyl file")
    argparser.add_argument("censat", type=str, help="Path to the CenSat BED file")
    argparser.add_argument("output", type=str, help="Output Path for the output files")

    # bed_parser arguments
    argparser.add_argument(
        "-m",
        "--mod_code",
        type=str,
        default="m",
        help='Modification code to filter bedMethyl file (default: "m")',
    )
    argparser.add_argument(
        "--methyl_bedgraph",
        action="store_true",
        default=False,
        help="Flag indicating if the input is a bedgraph. (default: False)",
    )
    argparser.add_argument(
        "--min_valid_cov",
        type=int,
        default=1,
        help="Minimum valid coverage to consider a methylation site (read from full modkit pileup files). (default: 10)",
    )
    argparser.add_argument(
        "-s",
        "--sat_type",
        type=str,
        default="active_hor",
        help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")',
    )
    argparser.add_argument(
        "--regions_prefiltered",
        action="store_true",
        default=False,
        help="Set flag if your annotations bed file is already subset to only the region you desire. (default: False)",
    )

    # calculate_matrices arguments
    argparser.add_argument(
        "--window_size",
        type=int,
        default=1190,
        help="Window size to calculate prior regions. (default: 1190)",
    )
    argparser.add_argument(
        "--step_size",
        type=int,
        default=1190,
        help="Step size when calculation windows for priors. (default: 1190)",
    )
    argparser.add_argument(
        "--prior_threshold",
        type=float,
        default=33.3,
        help="Threshold for determining if a window is a CDR. Uses this percentile if --prior_use_percentile is passed (default: 30.0)",
    )
    argparser.add_argument(
        "--prior_use_percentile",
        action="store_true",
        default=True,
        help="Whether or not to use percentile when calculating windowing priors. (default: True)",
    )
    argparser.add_argument(
        "--min_prior_size",
        type=int,
        default=11900,
        help="Minimum size for CDR regions. (default: 8330)",
    )
    argparser.add_argument(
        "--enrichment",
        action="store_true",
        default=False,
        help="Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)",
    )
    argparser.add_argument(
        "--percentile_emissions",
        action="store_true",
        default=False,
        help="Use values for flags w,x,y,z as raw threshold cutoffs for each emission category. (default: False)",
    )
    argparser.add_argument(
        "-x",
        type=float,
        default=25.0,
        help="Threshold of non-zero methylation percentile to be classified as low (default: 33.3)",
    )
    argparser.add_argument(
        "-y",
        type=float,
        default=50.0,
        help="Threshold of non-zero methylation percentile to be classified as medium (default: 66.6)",
    )
    argparser.add_argument(
        "-z",
        type=float,
        default=75.0,
        help="Threshold of non-zero methylation percentile to be classified as high (default: 100.0)",
    )

    # custom matrix input arguments
    argparser.add_argument(
        "--e_matrix",
        type=str,
        help="Custom Emission Matrix (Example: [[0.008,0.40,0.412,0.18],[0.025,0.9225,0.05,0.0025]]', 'chr1:[[0.008,0.40,0.412,0.18],[0.025,0.9225,0.05,0.0025]]...' or path to a file containing the matrix",
    )
    argparser.add_argument(
        "--t_matrix",
        type=str,
        help="Custom Transition Matrix (Example: '[[0.9999, 0.0001], [0.0025, 0.9975]]', 'chr1:[[0.9999, 0.0001], [0.0025, 0.9975]]...' or path to a file containing the matrix",
    )

    # HMM Flags
    argparser.add_argument(
        "--n_iter",
        type=int,
        default=1,
        help="Maximum number of iteration allowed for the HMM. (default: 1)",
    )
    argparser.add_argument(
        "--tol",
        type=float,
        default=10,
        help="Cutoff for model convergence in hmmlearn. (default: 10)",
    )
    argparser.add_argument(
        "--min_cdr_size",
        type=int,
        default=1190,
        help="Minimum size of region identified. (default: 1190)",
    )
    argparser.add_argument(
        "--min_cdr_score",
        type=float,
        default=90,
        help="The minimum HMM score [0-100] required to call a CDR. (default: 95)",
    )
    argparser.add_argument(
        "--min_low_conf_score",
        type=float,
        default=50,
        help="The minimum HMM score [0-100] required to call a low confidence CDR. (default: 75)",
    )
    argparser.add_argument(
        "--main_color",
        type=str,
        default="50,50,255",
        help="Color to dictate main regions. (default: 50,50,255)",
    )
    argparser.add_argument(
        "--low_conf_color",
        type=str,
        default="100,150,200",
        help="Color to dictate low confidence regions. (default: 100,150,200)",
    )

    # other_flags
    argparser.add_argument(
        "--min_subCDRs",
        type=int,
        default=3,
        help="Minimum number of subCDRs to report a CDR. (default: 3)",
    )
    argparser.add_argument(
        "--output_everything",
        action="store_true",
        default=False,
        help="Set to true if you would like to save all intermediate filesf. (default: False)",
    )
    argparser.add_argument(
        "--output_label",
        type=str,
        default="CDR",
        help='Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "subCDR")',
    )

    return argparser.parse_args()


def main():
    args = parse_command_line_arguments()

    # Validate matrix arguments
    if bool(args.e_matrix) != bool(args.t_matrix):
        raise ValueError("Both --e_matrix and --t_matrix must be provided together")

    output_prefix = os.path.splitext(args.output)[0]
    sat_types = [st.strip() for st in args.sat_type.split(",")]

    def generate_output_bed(all_chroms_dict, output_file, columns=["starts", "ends"]):
        all_lines = []
        for chrom in all_chroms_dict:
            chrom_data = all_chroms_dict[chrom]
            for i in range(len(chrom_data["starts"])):
                line = [chrom]
                for col in columns:
                    if col in chrom_data:
                        line.append(str(chrom_data[col][i])) 
                all_lines.append(line)
                
        all_lines = sorted(all_lines, key=lambda x: (x[0], int(x[1])))
        with open(output_file, 'w') as file:
            for line in all_lines: 
                file.write("\t".join(line) + "\n")

    # Parse in the bed files:
    parse = bed_parser(
        mod_code=args.mod_code,
        methyl_bedgraph=args.methyl_bedgraph,
        min_valid_cov=args.min_valid_cov,
        sat_type=sat_types,
        regions_prefiltered=args.regions_prefiltered,
    )

    regions_all_chroms, methylation_all_chroms = parse.process_files(
        methylation_path=args.bedmethyl, regions_path=args.censat
    )

    if args.output_everything:
        generate_output_bed(regions_all_chroms, f"{output_prefix}_regions.bed", columns=["starts", "ends"])
        generate_output_bed(methylation_all_chroms, f"{output_prefix}_methylation_in_regions.bedgraph", columns=["starts", "ends", "fraction_modified"])

    # create priors class as either way needs to be utilized
    priors = calculate_matrices(
        window_size=args.window_size,
        step_size=args.window_size,
        min_prior_size=args.min_prior_size,
        enrichment=args.enrichment,
        percentile_emissions=args.percentile_emissions,
        x=args.x,
        y=args.y,
        z=args.z,
        output_label=args.output_label,
    )

    (
        priors_all_chroms,
        windowmean_all_chroms,
        methylation_emissions_priors_all_chroms,
        emission_matrix_all_chroms,
        transition_matrix_all_chroms,
    ) = priors.priors_all_chromosomes(
        methylation_all_chroms=methylation_all_chroms,
        regions_all_chroms=regions_all_chroms,
        prior_threshold=args.prior_threshold,
    )

    if args.output_everything:
        generate_output_bed(priors_all_chroms, f"{output_prefix}_priors.bed", columns=["starts", "ends"])
        generate_output_bed(windowmean_all_chroms, f"{output_prefix}_windowmeans.bedgraph", columns=["starts", "ends", "means"])
        generate_output_bed(methylation_emissions_priors_all_chroms, f"{output_prefix}_emissions.bedgraph", columns=["starts", "ends", "emissions"])
        generate_output_bed(methylation_emissions_priors_all_chroms, f"{output_prefix}_priors.bedgraph", columns=["starts", "ends", "priors"])

        # Save matrices to a text file with key: numpy array representation
        with open(f"{output_prefix}_emission_matrices.txt", "w") as f:
            for key, matrix in emission_matrix_all_chroms.items():
                f.write(f"{key}: {matrix.tolist()}\n")

        with open(f"{output_prefix}_transition_matrices.txt", "w") as f:
            for key, matrix in transition_matrix_all_chroms.items():
                f.write(f"{key}: {matrix.tolist()}\n")

    CDRhmm = hmmCDR(
        n_iter=args.n_iter,
        tol=args.tol,
        min_cdr_size=args.min_cdr_size,
        min_cdr_score=args.min_cdr_score,
        min_low_conf_score=args.min_low_conf_score,
        main_color=args.main_color,
        low_conf_color=args.low_conf_color,
        output_label=args.output_label,
    )

    cdrs_all_chroms, methylation_emissions_priors_hmm_all_chroms = CDRhmm.hmm_all_chromosomes(
        methylation_emissions_priors_all_chroms=methylation_emissions_priors_all_chroms,
        emission_matrix_all_chroms=emission_matrix_all_chroms,
        transition_matrix_all_chroms=transition_matrix_all_chroms,
    )

    generate_output_bed(cdrs_all_chroms, f"{output_prefix}_hmmCDR.bed", columns=["starts","ends","name","score","strand","starts","ends","itemRgb"])
    if args.output_everything:
        generate_output_bed(methylation_emissions_priors_hmm_all_chroms, f"{output_prefix}_hmm_scores.bedgraph", columns=["starts","ends","hmm_cdr_score"])
        generate_output_bed(methylation_emissions_priors_hmm_all_chroms, f"{output_prefix}_cdr_assignments.bedgraph", columns=["starts","ends","cdr_calls"])

if __name__ == "__main__":
    main()
