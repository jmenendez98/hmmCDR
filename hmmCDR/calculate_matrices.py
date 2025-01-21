import argparse
import concurrent.futures
import os

import numpy as np
from scipy.stats import mannwhitneyu

from typing import Dict, Optional, List, Union

from hmmCDR.bed_parser import bed_parser


class calculate_matrices:
    def __init__(
        self,
        min_prior_size,
        enrichment,
        output_label,
        step_size=10 # how many CpGs on either side to include in prior calculation
    ):

        self.step_size = step_size
        self.min_prior_size = min_prior_size

        self.enrichment = enrichment
        self.output_label = output_label

    def calculate_regional_stats(self, methylation):
        methyl_starts = np.array(methylation["starts"], dtype=int)
        methyl_frac_mod = np.array(methylation["fraction_modified"], dtype=float)

        p_values = np.full(len(methyl_starts), 1.0, dtype=float)
        
        for i in range(len(methyl_starts)):
            start = max(0, i - self.step_size)
            end = min(len(methyl_starts), i + self.step_size + 1)

            current_region_frac_mods = methyl_frac_mod[start:end]

            _, p_value = mannwhitneyu(current_region_frac_mods, methyl_frac_mod, alternative="less", nan_policy="omit")
            p_values[i] = p_value

        methylation["mannU_p_value"] = p_values

        return methylation

    def find_priors(self, methylation, p_value_cutoff=0.01):
        methyl_starts = np.array(methylation["starts"], dtype=int)
        methyl_p_values = np.array(methylation["mannU_p_value"], dtype=float)
        methyl_sig = np.array([1 if p <= p_value_cutoff else 0 for p in methyl_p_values])

        # create the priors dictionary for the current chromosome
        priors = {"starts": [], "ends": []}

        priors_idx = np.where(methyl_sig == 1)[0]
        priors_idx_diff = np.diff(priors_idx)
        priors_idx_breaks = np.where(priors_idx_diff != 1)[0] + 1 
        priors_runs = np.split(priors_idx, priors_idx_breaks)

        for run in priors_runs:
            if len(run) == 0:
                continue

            min_idx, max_idx = np.min(run), np.max(run)
            start, end = methyl_starts[min_idx], methyl_starts[max_idx]+1
            run_length = end - start

            slice_scores = methyl_p_values[min_idx:max_idx+1]
            if slice_scores.size == 0:
                continue
            score = np.median(slice_scores)
            if np.isnan(score):
                continue

            if run_length > self.min_prior_size: 
                priors["starts"].append(start)
                priors["ends"].append(end)

        priors["starts"] = np.array(priors["starts"], dtype=int)
        priors["ends"] = np.array(priors["ends"], dtype=int)

        return priors

    def assign_priors(self, methylation, priors):
        priors_starts = np.array(priors['starts'], dtype=int)
        priors_ends = np.array(priors['ends'], dtype=int)
        methyl_starts = np.array(methylation['starts'], dtype=int)

        overlaps = np.zeros(len(methyl_starts), dtype=int)
        for region_start, region_end in zip(priors_starts, priors_ends):
            for i, methyl_start in enumerate(methyl_starts):
                if (methyl_start < region_end) and (methyl_start >= region_start):
                    overlaps[i] = 1

        if np.any(overlaps):
            methylation["priors"] = overlaps
        else: 
            ValueError(f"No priors found...")

        return methylation

    def assign_emissions(self, methylation):

        methyl_frac_mod = np.array(methylation["fraction_modified"], dtype=float)
        methylation_after_emissions_assigned = np.zeros(len(methyl_frac_mod), dtype=int)

        for i, score in enumerate(methyl_frac_mod):
            if score > self.x and score <= self.y:
                methylation_after_emissions_assigned[i] = 1
            elif score > self.y and score <= self.z:
                methylation_after_emissions_assigned[i] = 2
            elif score > self.z:
                methylation_after_emissions_assigned[i] = 3

        methylation["emissions"] = methylation_after_emissions_assigned

        return methylation

    def calculate_emission_matrix(self, methylation):
        # Create contingency table using NumPy for faster computation
        emission_matrix = np.zeros((2, 4))

        methyl_priors = np.array(methylation["priors"], dtype=int)
        methyl_emissions = np.array(methylation["emissions"], dtype=int)

        for state in [0, 1]:
            state_idx = np.where(methyl_priors == state)
            emissions_in_state = methyl_emissions[state_idx]

            for emission in emissions_in_state:
                emission_matrix[state, emission] += 1

            emission_matrix[state, :] = emission_matrix[state, :] / len(emissions_in_state)

        # Add a check to ensure rows sum to one (with a small tolerance for floating-point precision)
        assert np.allclose(emission_matrix.sum(axis=1), 1.0, rtol=1e-5), f"Transition matrix rows do not sum to one:\n{emission_matrix}"

        return emission_matrix

    def calculate_transition_matrix(self, methylation):
        # Convert priors to numpy array for faster operations
        priors = methylation["priors"] 

        # Create pairs of consecutive states using array slicing
        state_pairs = np.vstack((priors[:-1], priors[1:])).T

        # Use numpy's unique with return_counts to get transition counts
        unique_pairs, counts = np.unique(state_pairs, axis=0, return_counts=True)

        # Initialize transition matrix with zeros
        transition_matrix = np.zeros((2, 2))

        # Fill transition matrix using unique pairs and counts
        for (i, j), count in zip(unique_pairs, counts):
            transition_matrix[int(i), int(j)] = count

        # Normalize by row sums (avoiding division by zero)
        row_sums = transition_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        transition_matrix = transition_matrix / row_sums

        # Add a check to ensure rows sum to one (with a small tolerance for floating-point precision)
        assert np.allclose(transition_matrix.sum(axis=1), 1.0, rtol=1e-5), f"Transition matrix rows do not sum to one:\n{transition_matrix}"

        return transition_matrix

    def priors_single_chromosome(self, chrom, methylation, regions):
        methylation_mannu = self.calculate_regional_stats(methylation)

        # Find priors using current threshold
        priors = self.find_priors(methylation_mannu)

        # assign emissions to methylation bedtool
        methylation_mannu_emissions = self.assign_emissions(methylation)

        # add priors on to methylation bedtool
        methylation_mannu_emissions_priors = self.assign_priors(methylation=methylation_mannu_emissions, priors=priors)

        return ( chrom, priors, methylation_mannu_emissions_priors)

    def priors_all_chromosomes(self, methylation_all_chroms, regions_all_chroms):
        priors_all_chroms = {}
        methylation_emissions_priors_all_chroms = {}

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.priors_single_chromosome,
                    chrom, methylation_all_chroms[chrom], regions_all_chroms[chrom],
                ): chrom
                for chrom in methylation_all_chroms
            }

            for future in concurrent.futures.as_completed(futures):
                (
                    chrom,
                    priors,
                    methylation_emissions_priors,
                ) = future.result()

                priors_all_chroms[chrom] = priors
                methylation_emissions_priors_all_chroms[chrom] = methylation_emissions_priors

        return priors_all_chroms, methylation_emissions_priors_all_chroms


def main():
    argparser = argparse.ArgumentParser(
        description="Process bedMethyl and CenSat BED file to produce hmmCDR priors"
    )

    # required inputs
    argparser.add_argument("bedmethyl", type=str, help="Path to the bedmethyl file")
    argparser.add_argument("censat", type=str, help="Path to the censat BED file")
    argparser.add_argument("output", type=str, help="Path to the output BED file")

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
        help="Minimum valid coverage to consider a methylation site(read from full modkit pileup files). (default: 10)",
    )
    argparser.add_argument(
        "-s",
        "--sat_type",
        type=str,
        default="active_hor",
        help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "active_hor")',
    )
    argparser.add_argument(
        "--regions_prefiltered",
        action="store_true",
        default=False,
        help="Set flag if your annotations bed file is already subset to only the region you desire. (default: False)",
    )

    # calculate_matrices arguments
    argparser.add_argument(
        "--min_prior_size",
        type=int,
        default=1000,
        help="Minimum size for CDR regions. (default: 8330)",
    )
    argparser.add_argument(
        "--enrichment",
        action="store_true",
        default=False,
        help="Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)",
    )
    argparser.add_argument(
        "--output_label",
        type=str,
        default="subCDR",
        help='Label to use for name column of priorCDR BED file. (default: "subCDR")',
    )

    args = argparser.parse_args()
    sat_types = [st.strip() for st in args.sat_type.split(",")]
    output_prefix = os.path.splitext(args.output)[0]

    parse_beds = bed_parser(
        mod_code=args.mod_code,
        methyl_bedgraph=args.methyl_bedgraph,
        min_valid_cov=args.min_valid_cov,
        sat_type=sat_types,
        regions_prefiltered=args.regions_prefiltered,
    )

    regions_dict, methylation_dict = parse_beds.process_files(
        methylation_path=args.bedmethyl,
        regions_path=args.censat,
    )

    priors = calculate_matrices(
        min_prior_size=args.min_prior_size,
        enrichment=args.enrichment,
        output_label=args.output_label,
    )

    (
        priors_all_chroms,
        methylation_mannu_emissions_priors_all_chroms,
    ) = priors.priors_all_chromosomes(methylation_all_chroms=methylation_dict, regions_all_chroms=regions_dict)

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

    generate_output_bed(priors_all_chroms, f"{output_prefix}_priors.bed", columns=["starts", "ends"])
    generate_output_bed(methylation_mannu_emissions_priors_all_chroms, f"{output_prefix}_mannu.bedgraph", columns=["starts", "ends", "mannU_p_value"])

if __name__ == "__main__":
    main()