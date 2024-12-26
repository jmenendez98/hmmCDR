import argparse
import ast
import concurrent.futures
import os

import numpy as np
import pandas as pd
import pybedtools
from hmmlearn import hmm

from hmmCDR.bed_parser import bed_parser
from hmmCDR.calculate_matrices import calculate_matrices


class hmmCDR:
    def __init__(
        self,
        n_iter,
        tol,
        merge_distance,
        min_cdr_size,
        min_cdr_score,
        min_low_conf_size,
        min_low_conf_score,
        main_color,
        low_conf_color,
        output_label,
    ):

        self.n_iter = n_iter
        self.tol = tol
        self.merge_distance = merge_distance

        self.min_cdr_size = min_cdr_size
        self.min_cdr_score = min_cdr_score
        self.min_low_conf_size = min_low_conf_size
        self.min_low_conf_score = min_low_conf_score

        self.main_color = main_color
        self.low_conf_color = low_conf_color
        self.output_label = output_label

    def runHMM(self, labelled_methylation_data, transition_matrix, emission_matrix):
        """
        Apply a Hidden Markov Model (HMM) to analyze methylation data.

        This function uses a categorical HMM to decode methylation data and compute
        HMM scores based on emission probabilities. The resulting scores are added
        to the input DataFrame.

        Args:
            methylation_df (pd.DataFrame): A DataFrame containing methylation data with an 'emission' column.
            transition_matrix (np.ndarray): The transition probability matrix for the HMM.
            emission_matrix (np.ndarray): The emission probability matrix for the HMM.

        Returns:
            pd.DataFrame: The input DataFrame with an additional column, 'HMM_score', representing the computed HMM scores based on responsibilities.
        """
        model = hmm.CategoricalHMM(
            n_components=2, n_iter=self.n_iter, tol=self.tol, init_params=""
        )
        model.startprob_ = np.array([1.0, 0.0])
        model.transmat_ = transition_matrix
        model.emissionprob_ = emission_matrix

        emission_data = labelled_methylation_data["emission"].values.reshape(-1, 1)
        _, predicted_states = model.decode(emission_data, algorithm="viterbi")
        log_likelihood, responsibilities = model.score_samples(emission_data)

        labelled_methylation_data["HMM_score"] = [
            score[1] * 100 for score in responsibilities
        ]

        return labelled_methylation_data

    def create_subCDR_df(self, labelled_methylation_data):
        """
        Create a DataFrame of high-confidence and low-confidence CpG-dense regions (CDRs).

        Args:
            df (pd.DataFrame): Input DataFrame with columns 'chrom', 'start', 'end' as first three columns, and 'CDR_score' as a column in the dataframe.

        Returns:
            pd.DataFrame: Sorted DataFrame of annotated high-confidence and low-confidence CDRs.
            pd.DataFrame: Original DataFrame with scores used for creating CDRs.
        """
        # Step 1: Extract necessary columns
        hmmCDR_scores = labelled_methylation_data.iloc[:, :3]
        hmmCDR_scores["HMM_score"] = labelled_methylation_data[["HMM_score"]]

        def merge_CpG_dataframe(dataframe, threshold, merge_distance):
            """
            Merge CpG regions with scores above the threshold within a certain distance.

            Args:
                dataframe (pd.DataFrame): Input DataFrame with columns 'chrom', 'start', 'end', 'HMM_score'.
                threshold (float): Minimum score to include a CpG region.
                merge_distance (int): Maximum allowable gap between regions for merging.

            Returns:
                pd.DataFrame: Merged CpG regions with average scores.
            """
            # Filter rows by threshold
            filtered_df = dataframe[dataframe["HMM_score"] > threshold]
            if filtered_df.empty:
                return pd.DataFrame(columns=["chrom", "start", "end", "score"])

            # Convert to BedTool and merge nearby intervals
            bedtool = pybedtools.BedTool.from_dataframe(filtered_df)
            merged_bedtool = bedtool.merge(
                d=merge_distance, c=4, o="mean"
            )  # Use mean score

            # Convert back to DataFrame and round scores
            merged_df = merged_bedtool.to_dataframe(
                names=["chrom", "start", "end", "score"]
            )
            merged_df["score"] = merged_df["score"].round(5)

            return merged_df

        # Step 2: Generate high-confidence CDRs
        CDRs = merge_CpG_dataframe(
            hmmCDR_scores, self.min_cdr_score, self.merge_distance
        )

        # Filter regions by size
        CDRs["size"] = CDRs["end"] - CDRs["start"]
        CDRs = CDRs[CDRs["size"] >= self.min_cdr_size]
        # cdr_scores_column = CDRs['score']

        # Step 3: Handle low-confidence CDRs
        low_conf_CDRs = pd.DataFrame()
        if 0.0 < self.min_low_conf_score < self.min_cdr_score:
            low_conf_CDRs = merge_CpG_dataframe(
                hmmCDR_scores, self.min_low_conf_score, self.merge_distance
            )

            # Filter low-confidence regions by size
            low_conf_CDRs["size"] = low_conf_CDRs["end"] - low_conf_CDRs["start"]
            low_conf_CDRs = low_conf_CDRs[
                low_conf_CDRs["size"] >= self.min_low_conf_size
            ]

            # Subtract high-confidence regions from low-confidence regions
            cdr_bedtool = pybedtools.BedTool.from_dataframe(CDRs)
            low_conf_bedtool = pybedtools.BedTool.from_dataframe(low_conf_CDRs)
            low_conf_CDRs = low_conf_bedtool.subtract(cdr_bedtool).to_dataframe(
                names=["chrom", "start", "end", "score", "size"]
            )

            if not low_conf_CDRs.empty:
                low_conf_CDRs["name"] = f"low_conf_sub{self.output_label}"
                low_conf_CDRs["strand"] = "."

        # Step 4: Annotate high-confidence regions
        CDRs["name"] = f"sub{self.output_label}"
        CDRs["strand"] = "."

        # Combine and sort results
        try:
            combined_CDRs = pd.concat(
                [
                    CDRs[["chrom", "start", "end", "name", "score", "strand"]],
                    low_conf_CDRs[["chrom", "start", "end", "name", "score", "strand"]],
                ]
            )
        except KeyError:
            if low_conf_CDRs.empty:
                combined_CDRs = CDRs[
                    ["chrom", "start", "end", "name", "score", "strand"]
                ]
            elif CDRs.empty:
                combined_CDRs = low_conf_CDRs[
                    ["chrom", "start", "end", "name", "score", "strand"]
                ]
            else:
                combined_CDRs = pd.DataFrame(
                    columns=["chrom", "start", "end", "name", "score", "strand"]
                )

        combined_CDRs["thickStart"] = combined_CDRs["start"]
        combined_CDRs["thickEnd"] = combined_CDRs["end"]

        # Assign colors based on the 'name' column
        combined_CDRs["itemRgb"] = np.where(
            combined_CDRs["name"] == f"sub{self.output_label}",
            self.main_color,
            self.low_conf_color,
        )

        # Ensure integer types for relevant columns
        combined_CDRs = combined_CDRs.astype(
            {
                "start": "int64",
                "end": "int64",
                "thickStart": "int64",
                "thickEnd": "int64",
            }
        )

        # Sort by start position
        combined_CDRs = combined_CDRs.sort_values(by="start").reset_index(drop=True)

        return combined_CDRs, hmmCDR_scores

    def hmm_single_chromosome(
        self, chrom, labelled_methylation_data, emission_matrix, transition_matrix
    ):
        hmmlabelled_bed4Methyl = self.runHMM(
            labelled_methylation_data, transition_matrix, emission_matrix
        )
        hmmCDR_result, hmmCDR_scores = self.create_subCDR_df(hmmlabelled_bed4Methyl)
        return chrom, hmmCDR_result, hmmCDR_scores, emission_matrix, transition_matrix

    def hmm_all_chromosomes(
        self,
        labelled_methylation_chrom_dict,
        emission_matrix_chrom_dict=None,
        transition_matrix_chrom_dict=None,
    ):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.hmm_single_chromosome,
                    chrom,
                    labelled_methylation_chrom_dict[chrom],
                    emission_matrix_chrom_dict[chrom],
                    transition_matrix_chrom_dict[chrom],
                ): chrom
                for chrom in labelled_methylation_chrom_dict
            }

            results = {chrom: future.result() for future, chrom in futures.items()}
            hmmCDRresults_chrom_dict = {
                chrom: result[1] for chrom, result in results.items()
            }
            hmmCDRscores_chrom_dict = {
                chrom: result[2] for chrom, result in results.items()
            }

        return hmmCDRresults_chrom_dict, hmmCDRscores_chrom_dict


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
        "--bedgraph",
        action="store_true",
        default=False,
        help="Flag indicating if the input is a bedgraph. (default: False)",
    )
    argparser.add_argument(
        "--min_valid_cov",
        type=int,
        default=10,
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
        "--pre_subset_censat",
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
        default=20.0,
        help="Threshold for determining if a window is a CDR. Uses this percentile if --prior_use_percentile is passed (default: 30.0)",
    )
    argparser.add_argument(
        "--prior_use_percentile",
        action="store_true",
        default=False,
        help="Whether or not to use percentile when calculating windowing priors. (default: False)",
    )
    argparser.add_argument(
        "--min_prior_size",
        type=int,
        default=8330,
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
        "-w",
        type=float,
        default=0.0,
        help="Threshold of non-zero methylation percentile to be classified as None (default: 0.0)",
    )
    argparser.add_argument(
        "-x",
        type=float,
        default=33.3,
        help="Threshold of non-zero methylation percentile to be classified as low (default: 33.3)",
    )
    argparser.add_argument(
        "-y",
        type=float,
        default=66.6,
        help="Threshold of non-zero methylation percentile to be classified as medium (default: 66.6)",
    )
    argparser.add_argument(
        "-z",
        type=float,
        default=90.0,
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
        "--hmm_merge_distance",
        type=int,
        default=1190,
        help="Distance to merge adjacently labelled subCDR regions. (default: 1190)",
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
        default=99,
        help="The minimum HMM score [0-100] required to call a CDR. (default: 95)",
    )
    argparser.add_argument(
        "--min_low_conf_size",
        type=int,
        default=0,
        help="Minimum size of region identified. (default: 0)",
    )
    argparser.add_argument(
        "--min_low_conf_score",
        type=float,
        default=75,
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
        "--large_merge_distance",
        type=int,
        default=200000,
        help="Distance to merge subCDRs into a larger CDR annotation. (default: 200000)",
    )
    argparser.add_argument(
        "--output_all",
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

    # Parse in the bed files:
    parseCDRBeds = bed_parser(
        mod_code=args.mod_code,
        bedgraph=args.bedgraph,
        min_valid_cov=args.min_valid_cov,
        sat_type=sat_types,
        pre_subset_censat=args.pre_subset_censat,
    )

    methylation_chrom_dict, regions_chrom_dict = parseCDRBeds.process_files(
        bedmethyl_path=args.bedmethyl, censat_path=args.censat
    )

    # create priors class as either way needs to be utilized
    priors = calculate_matrices(
        window_size=args.window_size,
        step_size=args.window_size,
        min_prior_size=args.min_prior_size,
        enrichment=args.enrichment,
        percentile_emissions=args.percentile_emissions,
        w=args.w,
        x=args.x,
        y=args.y,
        z=args.z,
        output_label=args.output_label,
    )
    labelled_methylation_chrom_dict = {}

    # Skip prior finding if matrices are provided
    if args.e_matrix and args.t_matrix:
        # Parse emission matrix from string or file
        emission_matrix_chrom_dict = {}
        try:
            if os.path.isfile(args.e_matrix):
                with open(args.e_matrix, "r") as f:
                    emission_input = f.read().strip()
            else:
                emission_input = args.e_matrix

            # Check for label
            if ":" in emission_input:
                chrom, matrix_str = emission_input.split(":", 1)
                emission_matrix_chrom_dict[chrom.strip()] = ast.literal_eval(
                    matrix_str.strip()
                )
            else:
                # No chromosome specified, use for all chromosomes in methylation_chrom_dict
                matrix = ast.literal_eval(emission_input)
                for chrom in methylation_chrom_dict:
                    emission_matrix_chrom_dict[chrom] = matrix

        except (ValueError, SyntaxError, FileNotFoundError):
            # [[0.008791358571574608, 0.4036358060449109, 0.4158005929055781, 0.1717722424779364], [0.025597269624573378, 0.9215017064846417, 0.05005688282138794, 0.002844141069397042]]
            raise ValueError(
                "Invalid emission matrix format. Use format like: '[[0.008,0.40,0.412,0.18],[0.025,0.9225,0.05,0.0025]]', 'chr1:[[0.008,0.40,0.412,0.18],[0.025,0.9225,0.05,0.0025]]...' or path to a file containing the matrix"
            )

        # Check for missing chromosome matrices and print warnings
        for chrom in methylation_chrom_dict:
            if chrom not in emission_matrix_chrom_dict:
                print(
                    f"Warning: No emission matrix given for chromosome {chrom}.\n\tProceeding with default emission matrix for {chrom}: [[0.008,0.40,0.412,0.18],[0.025,0.9225,0.05,0.0025]] (not recommended if emissions are changed!!!)"
                )
                emission_matrix_chrom_dict[chrom] = [
                    [0.008, 0.40, 0.412, 0.18],
                    [0.025, 0.9225, 0.05, 0.0025],
                ]

        # Parse transition matrix from string or file
        transition_matrix_chrom_dict = {}
        try:
            if os.path.isfile(args.t_matrix):
                with open(args.t_matrix, "r") as f:
                    transition_input = f.read().strip()
            else:
                transition_input = args.t_matrix

            # Check for label
            if ":" in transition_input:
                chrom, matrix_str = transition_input.split(":", 1)
                transition_matrix_chrom_dict[chrom.strip()] = ast.literal_eval(
                    matrix_str.strip()
                )
            else:
                # No chromosome specified, use for all chromosomes in methylation_chrom_dict
                matrix = ast.literal_eval(transition_input)
                for chrom in methylation_chrom_dict:
                    transition_matrix_chrom_dict[chrom] = matrix

        except (ValueError, SyntaxError, FileNotFoundError):
            raise ValueError(
                "Invalid transition matrix format. Use format like: '[[0.9999, 0.0001], [0.0025, 0.9975]]', 'chr1:[[0.9999, 0.0001], [0.0025, 0.9975]]...' or path to a file containing the matrix"
            )

        # Check for missing chromosome transition matrices and print warnings
        for chrom in methylation_chrom_dict:
            if chrom not in transition_matrix_chrom_dict:
                print(
                    f"Warning: No transition matrix given for chromosome {chrom}.\n\tProceeding with default transition matrix for {chrom}: [[0.9999, 0.0001], [0.0025, 0.9975]]"
                )
                transition_matrix_chrom_dict[chrom] = [
                    [0.9999, 0.0001],
                    [0.0025, 0.9975],
                ]

        # label methylation data if matrices are passed in...
        for chrom in methylation_chrom_dict:
            labelled_methylation_chrom_dict[chrom] = priors.assign_emissions(
                methylation_chrom_dict[chrom],
                priors.calculate_emission_thresholds(methylation_chrom_dict[chrom]),
            )

    else:
        (
            priors_chrom_dict,
            windowmean_chrom_dict,
            labelled_methylation_chrom_dict,
            emission_matrix_chrom_dict,
            transition_matrix_chrom_dict,
        ) = priors.priors_all_chromosomes(
            methylation_chrom_dict=methylation_chrom_dict,
            regions_chrom_dict=regions_chrom_dict,
            prior_percentile=args.prior_use_percentile,
            prior_threshold=args.prior_threshold,
        )

        if args.output_all:
            # Concatenate, sort and save bed files
            combined_priors = pd.concat(priors_chrom_dict.values(), ignore_index=True)
            combined_priors = combined_priors.sort_values(
                by=combined_priors.columns[:2].tolist()
            )
            combined_priors.to_csv(
                f"{output_prefix}_priors.bed", sep="\t", index=False, header=False
            )

            combined_windowmean = pd.concat(
                windowmean_chrom_dict.values(), ignore_index=True
            )
            combined_windowmean = combined_windowmean.sort_values(
                by=combined_windowmean.columns[:2].tolist()
            )
            combined_windowmean.to_csv(
                f"{output_prefix}_windowmean.bedgraph",
                sep="\t",
                index=False,
                header=False,
            )

            # Save matrices to a text file with key: numpy array representation
            with open(f"{output_prefix}_emission_matrices.txt", "w") as f:
                for key, matrix in emission_matrix_chrom_dict.items():
                    f.write(f"{key}: {matrix.tolist()}\n")

            with open(f"{output_prefix}_transition_matrices.txt", "w") as f:
                for key, matrix in transition_matrix_chrom_dict.items():
                    f.write(f"{key}: {matrix.tolist()}\n")

    CDRhmm = hmmCDR(
        n_iter=args.n_iter,
        tol=args.tol,
        merge_distance=args.hmm_merge_distance,
        min_cdr_size=args.min_cdr_size,
        min_cdr_score=args.min_cdr_score,
        min_low_conf_size=args.min_low_conf_size,
        min_low_conf_score=args.min_low_conf_score,
        main_color=args.main_color,
        low_conf_color=args.low_conf_color,
        output_label=args.output_label,
    )

    hmm_results_chrom_dict, hmm_scores_chrom_dict = CDRhmm.hmm_all_chromosomes(
        labelled_methylation_chrom_dict=labelled_methylation_chrom_dict,
        emission_matrix_chrom_dict=emission_matrix_chrom_dict,
        transition_matrix_chrom_dict=transition_matrix_chrom_dict,
    )

    # output subCDR
    combined_subCDRs = pd.concat(hmm_results_chrom_dict.values(), axis=0)
    combined_subCDRs = combined_subCDRs.sort_values(
        by=list(combined_subCDRs.columns[:2])
    )
    combined_subCDRs.to_csv(
        f"{output_prefix}_sub{args.output_label}.bed",
        sep="\t",
        index=False,
        header=False,
    )

    if args.output_all:  # output HMM score bedgraph
        combined_hmm_scores = pd.concat(hmm_scores_chrom_dict.values(), axis=0)
        combined_hmm_scores = combined_hmm_scores.sort_values(
            by=list(combined_hmm_scores.columns[:2])
        )
        combined_hmm_scores.to_csv(
            f"{output_prefix}_scores.bedgraph", sep="\t", index=False, header=False
        )

    # create final CDR output that merges adjacent subCDRs and reports scoring of over 3
    combined_subCDR_bedtool = pybedtools.BedTool.from_dataframe(
        combined_subCDRs[combined_subCDRs.iloc[:, 3] == f"sub{args.output_label}"]
    )
    combined_big_CDR = combined_subCDR_bedtool.merge(
        d=args.large_merge_distance, c=2, o="count"
    ).to_dataframe(names=["chrom", "start", "end", "count"])
    combined_big_CDR["name"] = combined_big_CDR.loc[:, "count"].apply(
        lambda x: (
            f"low_conf_{args.output_label}"
            if x < args.min_subCDRs
            else f"{args.output_label}"
        )
    )

    # Format final output dataframe
    bigCDR_names = combined_big_CDR.loc[:, "name"]
    bigCDR_counts = combined_big_CDR.loc[:, "count"]
    combined_big_CDR = combined_big_CDR.iloc[:, :3]
    combined_big_CDR["name"] = bigCDR_names
    combined_big_CDR["score"] = bigCDR_counts
    combined_big_CDR["strand"] = "."
    combined_big_CDR["thickStart"] = combined_big_CDR.iloc[:, 1]
    combined_big_CDR["thickEnd"] = combined_big_CDR.iloc[:, 2]
    combined_big_CDR["itemRgb"] = np.where(
        combined_big_CDR["name"] == f"{args.output_label}",
        args.main_color,
        np.where(
            combined_big_CDR["name"] == f"low_conf_{args.output_label}",
            args.low_conf_color,
            "",
        ),
    )

    # Sort the final output
    combined_big_CDR = combined_big_CDR.sort_values(by=["chrom", "start"])
    combined_big_CDR.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
