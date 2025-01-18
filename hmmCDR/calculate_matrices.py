import argparse
import concurrent.futures
import os

import numpy as np
import pandas as pd

from typing import Dict, Optional, List, Union

from hmmCDR.bed_parser import bed_parser


class calculate_matrices:
    def __init__(
        self,
        window_size,
        step_size,
        min_prior_size,
        percentile_emissions,
        enrichment,
        output_label,
        x=25, # anything under 25% methylation
        y=50, # anything between 25% and 75% methylation
        z=75, # anything over 75% methylation
    ):

        self.window_size = window_size
        self.step_size = step_size
        self.min_prior_size = min_prior_size

        self.percentile_emissions = percentile_emissions

        self.enrichment = enrichment
        self.output_label = output_label

        self.x, self.y, self.z = x, y, z

    def create_windows(self, regions):
        windows: Dict[str, list] = {}

        windows["starts"], windows["ends"] = [], []

        for start, end in zip(regions["starts"], regions["ends"]):
            starts = np.arange(start, end-self.window_size, self.step_size, dtype=int)
            ends = starts + self.window_size

            windows["starts"].extend(starts)
            windows["ends"].extend(ends)

        return windows

    def mean_methylation_in_windows(self, methylation, windows):

        windows["means"] = np.empty(len(windows["starts"]), dtype=float)
        methyl_starts = np.array(methylation["starts"], dtype=int)
        methyl_scores = np.array(methylation["scores"], dtype=float)

        for i, (region_start, region_end) in enumerate(zip(windows["starts"], windows["ends"])):
            overlaps = np.where(np.logical_and(methyl_starts<region_end, methyl_starts>=region_start))

            if len(methyl_scores[overlaps]) > 0:
                # Calculate mean of scores for overlapping sites
                windows["means"][i] = np.mean(methyl_scores[overlaps])
            else:
                windows["means"][i] = np.nan

        return windows

    def get_prior_threshold(self, window_means, percentile):
        window_means_no_na = window_means["means"].dropna()
        return np.percentile(window_means_no_na, q=percentile)

    def find_priors(self, windows_means, threshold):
        prior_windows = [(windows_means["starts"][i], windows_means["ends"][i]) for i, window in enumerate(windows_means["means"]) if window < threshold]
        
        prior_windows = sorted(prior_windows)
        merged = []

        current_start, current_end = prior_windows[0]
        for start, end in prior_windows[1:]:
            # Check if current interval overlaps or is adjacent to the next
            if start <= current_end + 1:
                # Merge by extending current_end if necessary
                current_end = max(current_end, end)
            else:
                # No overlap or adjacency, add current interval and start new one
                merged.append((current_start, current_end))
                current_start, current_end = start, end

        # Add the last interval
        merged.append((current_start, current_end))

        # filter out priors smaller than size threshold
        merged_filtered = [(start, end) for start, end in merged if end - start >= self.min_prior_size]

        # create the priors dictionary for the current chromosome
        priors = {"starts": [], "ends": []}
        for start, end in merged_filtered:
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

        methylation_scores = np.array(methylation["scores"], dtype=float)
        methylation_after_emissions_assigned = np.zeros(len(methylation_scores), dtype=int)

        for i, score in enumerate(methylation_scores):
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

    def priors_single_chromosome(self, chrom, methylation, regions, prior_threshold):
        windows = self.create_windows(regions)
        windows_means = self.mean_methylation_in_windows(methylation, windows)

        # Find priors using current threshold
        priors = self.find_priors(windows_means, prior_threshold)

        # assign emissions to methylation bedtool
        methylation_emissions = self.assign_emissions(methylation)

        if all(emission == 0 for emission in methylation_emissions["emissions"]):
            # If no priors found, adjust threshold
            print(f"No prior CDRs Found on {chrom} with prior threshold - {prior_threshold}")
            print(f"Continuing with defaults:\n emission matrix: [[0.002,0.10,0.28,0.60],[0.05,0.85,0.08,0.02]]\ntransition matrix: [[0.9999,0.003],[0.0001,0.997]]")
            
            emission_matrix = np.array([[0.008, 0.40, 0.412, 0.18], [0.025, 0.9225, 0.05, 0.0025]])
            transition_matrix = np.array([[0.9999, 0.0001], [0.0025, 0.9975]])

            return (
                chrom,
                priors,
                windows_means,
                methylation_emissions,
                emission_matrix,
                transition_matrix,
            )

        # add priors on to methylation bedtool
        methylation_emissions_priors = self.assign_priors(methylation=methylation_emissions, priors=priors)

        # calculate emission and transition matrices with assigned priors
        emission_matrix = self.calculate_emission_matrix(methylation_emissions_priors)
        transition_matrix = self.calculate_transition_matrix(methylation_emissions_priors)

        return ( chrom, priors, windows_means, methylation_emissions_priors, emission_matrix, transition_matrix )

    def priors_all_chromosomes(self, methylation_all_chroms, regions_all_chroms, prior_threshold):
        priors_all_chroms = {}
        windowmeans_all_chroms = {}
        methylation_emissions_priors_all_chroms = {}
        emission_matrix_all_chroms = {}
        transition_matrix_all_chroms = {}
        chromosomes = methylation_all_chroms.keys()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.priors_single_chromosome,
                    chrom, methylation_all_chroms[chrom], regions_all_chroms[chrom], prior_threshold,
                ): chrom
                for chrom in chromosomes
            }

            for future in concurrent.futures.as_completed(futures):
                (
                    chrom,
                    priors,
                    window_means,
                    methylation_emissions_priors,
                    emission_matrix,
                    transition_matrix,
                ) = future.result()

                priors_all_chroms[chrom] = priors
                windowmeans_all_chroms[chrom] = window_means
                methylation_emissions_priors_all_chroms[chrom] = methylation_emissions_priors
                emission_matrix_all_chroms[chrom] = emission_matrix
                transition_matrix_all_chroms[chrom] = transition_matrix

        return priors_all_chroms, windowmeans_all_chroms, methylation_emissions_priors_all_chroms, emission_matrix_all_chroms, transition_matrix_all_chroms


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
        default=10,
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
        help="Threshold for determining if a window is a CDR. Uses value as a percentile of all windows if --prior_use_percentile is passed (default: 30.0)",
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
        "--output_label",
        type=str,
        default="subCDR",
        help='Label to use for name column of priorCDR BED file. (default: "subCDR")',
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
        default=25,
        help="Threshold of non-zero methylation percentile to be classified as low (default: 33.3)",
    )
    argparser.add_argument(
        "-y",
        type=float,
        default=50,
        help="Threshold of non-zero methylation percentile to be classified as medium (default: 66.6)",
    )
    argparser.add_argument(
        "-z",
        type=float,
        default=75,
        help="Threshold of non-zero methylation percentile to be classified as high (default: 100.0)",
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
        window_size=args.window_size, step_size=args.window_size,
        min_prior_size=args.min_prior_size, percentile_emissions=args.percentile_emissions,
        enrichment=args.enrichment,
        x=args.x, y=args.y, z=args.z,
        output_label=args.output_label,
    )

    (
        priors_all_chroms,
        windowmean_all_chroms,
        methylation_emissions_priors_all_chroms,
        emission_matrix_all_chroms,
        transition_matrix_all_chroms,
    ) = priors.priors_all_chromosomes(methylation_all_chroms=methylation_dict, regions_all_chroms=regions_dict, prior_threshold=args.prior_threshold)

    def generate_output_bed(all_chroms_dict, output_file, columns=["starts", "ends"]):
        with open(output_file, 'w') as file: 
            for chrom in all_chroms_dict:
                chrom_data = all_chroms_dict[chrom]
                for i in range(len(all_chroms_dict["starts"])):
                    for col in columns:
                       if col in chrom_data:
                            line.append(str(chrom_data[col][i])) 
                file.write("\t".join(line) + "\n")

    generate_output_bed(priors_all_chroms, f"{args.output_prefix}_priors.bed", columns=["starts", "ends"])
    generate_output_bed(priors_all_chroms, f"{args.output_prefix}_windowmean.bedgraph", columns=["starts", "ends"])
    
    # Save matrices to a text file with key: numpy array representation
    with open(f"{output_prefix}_emission_matrices.txt", "w") as file:
        for key, matrix in emission_matrix_all_chroms.items():
            file.write(f"{key}: {matrix.tolist()}\n")

    with open(f"{output_prefix}_transition_matrices.txt", "w") as file:
        for key, matrix in transition_matrix_all_chroms.items():
            file.write(f"{key}: {matrix.tolist()}\n")


if __name__ == "__main__":
    main()
