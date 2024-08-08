import pandas as pd
import numpy as np 
import pybedtools
import argparse
from hmmlearn import hmm
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Process input files with optional parameters.')

    # Required arguments
    parser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    parser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    parser.add_argument('output_path', type=str, help='Path to the output priorCDRs BED file')

    # Create mutually exclusive group for file input options
    inputs_group = parser.add_mutually_exclusive_group(required=True)
    # Option for providing a single BED file
    inputs_group.add_argument('priors', defult=None,type=str, help='Path to the priorCDR bedfile')
    # Option for providing two matrix TSV files
    inputs_group.add_argument('-e', '--emission_matrix', default=None, type=str, help='Path to the emission matrix TSV file')
    inputs_group.add_argument('-t', '--transition_matrix', default=None, type=str, help='Path to the transition matrix TSV file')

    # Optional flags with default values
    parser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    parser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Satellite type/name to filter CenSat bed file. (default: "H1L")')
    parser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "CDR")')
    parser.add_argument('--use_percentiles', action='store_true', default=True, help='Use values for flags w,x,y,z as percentile cutoffs for each category. (default: False)')
    parser.add_argument('--n_iter', type=int, default=1, help='Maximum number of iteration allowed for the HMM. (default: 1)')
    parser.add_argument('-w', type=int, default=0, help='Theshold for methylation to be classified as very low (default: 0)')
    parser.add_argument('-x', type=int, default=25, help='Theshold for methylation to be classified as low (default: 25)')
    parser.add_argument('-y', type=int, default=50, help='Theshold for methylation to be classified as medium (default: 50)')
    parser.add_argument('-z', type=int, default=75, help='Theshold for methylation to be classified as high (default: 75)')

    args = parser.parse_args()
    return args

class hmmCDR:
    def __init__(self, bedMethyl_path, cenSat_path, output_path,
                 mod_code, sat_type, output_label, n_iter,
                 priors, emission_matrix, transition_matrix,
                 use_percentiles, w, x, y, z):

        self.mod_code = mod_code
        self.sat_type = sat_type
        self.output_label = output_label
        self.n_iter = n_iter
        self.output_path = output_path

        bedgraphMethyl = self.filter_bedMethyl(bedMethyl_path)
        filtered_regions = self.filter_regions(cenSat_path)
        intersected_bedgraphMethyl = self.intersect_files(bedgraphMethyl, filtered_regions)
            
        self.use_percentiles = use_percentiles
        if self.use_percentiles is True:
            self.w = self.calculate_percentiles(intersected_bedgraphMethyl, w)
            self.x = self.calculate_percentiles(intersected_bedgraphMethyl, x)
            self.y = self.calculate_percentiles(intersected_bedgraphMethyl, y)
            self.z = self.calculate_percentiles(intersected_bedgraphMethyl, z)
        else:
            self.w = w
            self.x = x
            self.y = y
            self.z = z
        
        if priors is not None:
            self.priors_df = pd.read_csv(priors, sep='\t', header=None)
            prior_labelled_bedgraphMethyl = self.assign_priors(intersected_bedgraphMethyl)
            emission_prior_labelled_bedgraphMethyl = self.assign_emissions(prior_labelled_bedgraphMethyl)

            self.emission_matrix = self.calculate_emission_matrix(emission_prior_labelled_bedgraphMethyl)
            self.transition_matrix = self.calculate_transition_matrix(emission_prior_labelled_bedgraphMethyl)
        elif emission_matrix is not None and transition_matrix is not None:
            self.emission_matrix = emission_matrix
            self.transition_matrix = transition_matrix
            self.priors = None
        else:
            raise ValueError("Invalid input: Provide either a priors file or both emission and transition matrix files.")
        
        emission_labelled_bedgraphMethyl = self.assign_emissions(intersected_bedgraphMethyl)
        predicted_state_labelled_bedgraphMethyl = self.runHMM(emission_labelled_bedgraphMethyl)
        
        hmmCDR_df = self.create_hmmCDRbed(predicted_state_labelled_bedgraphMethyl)

        
    def filter_bedMethyl(self, bedMethyl_path):
        df = pd.read_csv(bedMethyl_path, sep='\t', header=None, usecols=[0, 1, 2, 3, 10])
        filtered_df = df[df[3] == self.mod_code]
        filtered_df = filtered_df.drop(columns=[3])
        return filtered_df

    def filter_regions(self, cenSat_path):
        regions_df = pd.read_csv(cenSat_path, sep='\t', header=None)
        filtered_regions_df = regions_df[regions_df[3].str.contains(self.sat_type)]
        return filtered_regions_df

    def intersect_files(self, bedgraphMethyl, filtered_regions):
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bedgraphMethyl)
        regions_bedtool = pybedtools.BedTool.from_dataframe(filtered_regions)
        intersected = bedMethyl_bedtool.intersect(regions_bedtool, wa=True, u=True)
        intersected_df = intersected.to_dataframe(names=[0, 1, 2, 3])
        return intersected_df
    
    def calculate_percentiles(self, intersected_bedgraphMethyl, percentile):
        intersected_bedgraphMethyl[3].replace('.', np.nan, inplace=True)
        methylation_scores = intersected_bedgraphMethyl[3].dropna().astype(float)
        value = np.percentile(methylation_scores, q=self.percentile)
        return value
    
    def assign_priors(self, intersected_bedgraphMethyl):
        # Convert DataFrames to BedTool objects
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(intersected_bedgraphMethyl)
        priorCDR_df = self.priors_df[self.priors_df[3] == f'{self.output_label}'].drop(columns=[3])
        priorCDR_bedtool = pybedtools.BedTool.from_dataframe(priorCDR_df)
        priorTransition_df = self.priors_df[self.priors_df[3] == f'{self.output_label}_transition'].drop(columns=[3])
        priorTransition_bedtool = pybedtools.BedTool.from_dataframe(priorTransition_df)

        # Create a column for CDR state, default to 'A' (Not in CDR)
        intersected_bedgraphMethyl['state'] = 'A'

        # Label CDR regions as 'C'
        intersected_cdr = bedMethyl_bedtool.intersect(priorCDR_bedtool, wa=True, wb=True)
        cdr_df = intersected_cdr.to_dataframe()
        # Update the state based on matching 'start' values
        intersected_bedgraphMethyl.loc[intersected_bedgraphMethyl['start'].isin(cdr_df['start']), 'state'] = 'C'

        # Label CDR transition regions as 'B'
        intersected_transition = bedMethyl_bedtool.intersect(priorTransition_bedtool, wa=True, wb=True)
        transition_df = intersected_transition.to_dataframe()
        # Update the state based on matching 'start' values
        intersected_bedgraphMethyl.loc[intersected_bedgraphMethyl['start'].isin(transition_df['start']), 'state'] = 'B'

        return intersected_bedgraphMethyl

    def assign_priors(self, intersected_bedgraphMethyl):
        def assign_emission(value):
            if value > self.z:
                return 3
            elif value > self.y:
                return 2
            elif value > self.x:
                return 1
            elif value >= self.w:
                return 0
            
        intersected_bedgraphMethyl['emission'] = intersected_bedgraphMethyl['name'].apply(assign_emission)
        
        return intersected_bedgraphMethyl

    def calculate_transition_matrix(self, labeled_bedMethyl_df):
        transitions = {
            'A->A': 0,
            'A->B': 0,
            'A->C': 0,
            'B->A': 0,
            'B->B': 0,
            'B->C': 0,
            'C->A': 0,
            'C->B': 0,
            'C->C': 0
        }

        previous_state = labeled_bedMethyl_df.iloc[0]['state']
        
        for i in range(1, len(labeled_bedMethyl_df)):
            current_state = labeled_bedMethyl_df.iloc[i]['state']
            
            transitions[f'{previous_state}->{current_state}'] += 1
                
            previous_state = current_state
        
        total_A = transitions['A->A'] + transitions['A->B'] + transitions['A->C']
        total_B = transitions['B->A'] + transitions['B->B'] + transitions['B->C']
        total_C = transitions['C->A'] + transitions['C->B'] + transitions['C->C']
        
        transition_matrix = [
            [transitions['A->A'] / total_A if total_A else 0, transitions['A->B'] / total_A if total_A else 0, transitions['A->C'] / total_A if total_A else 0],
            [transitions['B->A'] / total_B if total_B else 0, transitions['B->B'] / total_B if total_B else 0, transitions['B->C'] / total_B if total_B else 0],
            [transitions['C->A'] / total_C if total_C else 0, transitions['C->B'] / total_C if total_C else 0, transitions['C->C'] / total_C if total_C else 0]
        ]
        
        return transition_matrix

    def calculate_emission_matrix(self, labeled_bedMethyl):
        # Get the counts of each emission in each state
        emission_counts = labeled_bedMethyl.groupby(['state', 'emission']).size().unstack(fill_value=0)
        
        # Normalize counts to probabilities to get the emission matrix
        emission_matrix = emission_counts.div(emission_counts.sum(axis=1), axis=0)
        
        return emission_matrix

    def runHMM(self, emission_labelled_bedgraphMethyl):
        """
        Create an HMM model with specified emission and transition matrices. Then Fit the HMM model on the emission data.

        Parameters:
        emission_matrix (numpy.ndarray): The emission probabilities matrix.
        transition_matrix (numpy.ndarray): The transition probabilities matrix.
        n_states (int): Number of states in the HMM.

        Returns:

        numpy.ndarray: The predicted states.
        """
        # Create an HMM instance with Multinomial emissions
        model = hmm.CategoricalHMM(n_components=self.n_states, n_features=4, n_iter=self.n_iter, init_params="")
        # Set the start probs
        model.startprob_ = np.array([1.0, 0.0, 0.0])
        # Set the transition matrix
        model.transmat_ = self.transition_matrix
        print('TransMat:', model.transmat_.shape, '\n', model.transmat_)
        # Set the emission matrix
        model.emissionprob_ = self.emission_matrix
        print('EmisMat:', model.emissionprob_.shape, '\n', model.emissionprob_)
        emission_data = emission_labelled_bedgraphMethyl['emission'].to_numpy().reshape(-1, 1)
        # Predict the hidden states
        logprob, predicted_states = model.decode(emission_data, algorithm="viterbi")
        emission_labelled_bedgraphMethyl['predicted_state'] = predicted_states
        return emission_labelled_bedgraphMethyl

    def create_hmmCDR_df(self, df):
        def merge_and_label_regions(df, state_value, label, merge_distance=1000):
            # Extract regions with the specified state
            state_df = df[df['predicted_state'] == state_value].copy()
            # Ensure columns are correctly named for pybedtools
            state_df.columns = ['chrom', 'start', 'end', 'name', 'state', 'emission', 'predicted_state']
            # Convert to a BedTool object
            state_bedtool = pybedtools.BedTool.from_dataframe(state_df[['chrom', 'start', 'end', 'name']])
            # Merge adjacent entries within the specified distance
            merged_bedtool = state_bedtool.merge(d=merge_distance)
            # Convert back to DataFrame and add the label column
            merged_df = merged_bedtool.to_dataframe(names=['chrom', 'start', 'end', 'name'])
            merged_df['name'] = label
            return merged_df
        
        merged_cdr_df = merge_and_label_regions(df, 2, f'{self.output_label}', 1000)
        merged_transition_df = merge_and_label_regions(df, 1, f'{self.output_label}_transition', 1000)

        # Combine the two DataFrames
        combined_df = pd.concat([merged_cdr_df, merged_transition_df], ignore_index=True)
        combined_df = combined_df.sort_values(by=['chrom', 'start']).reset_index(drop=True)

        def fix_transitions(df, distance=1000):
            df = df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
            for i in range(len(df)):
                if df.iloc[i]['name'] == 'CDR_transition':
                    prev_row = df.iloc[i - 1] if i > 0 else None
                    next_row = df.iloc[i + 1] if i < len(df) - 1 else None
                    if prev_row is not None and prev_row['name'] == 'CDR':
                        # Check distance from transition start to CDR end
                        if df.iloc[i]['start'] <= prev_row['end'] + distance:
                            # Adjust start and end to be adjacent
                            df.at[i, 'start'] = prev_row['end'] + 1
                    if next_row is not None and next_row['name'] == 'CDR':
                        # Check distance from transition end to CDR start
                        if next_row['start'] <= df.iloc[i]['end'] + distance:
                            # Adjust start and end to be adjacent
                            df.at[i, 'end'] = next_row['start'] - 1
            # Remove any transitions that end before they start (possible due to adjustment)
            df = df[df['start'] <= df['end']]
            return df
        
        adjusted_df = fix_transitions(combined_df, 1000)

        return adjusted_df


def main():
    args = parse_args()

    hmmCDR_generator = hmmCDR(args.bedMethyl_path,
                              args.cenSat_path,
                              args.output_path,
                              args.mod_code,
                              args.sat_type,
                              args.output_label,
                              args.n_iter,
                              args.priors,
                              args.emission_matrix,
                              args.transition_matrix,
                              args.use_percentiles,
                              args.w, args.x, args.y, args.z)

    hmmCDR_generator.hmmCDR_df.to_csv(f'{args.output_path}', sep='\t', header=False, index=False)

if __name__ == "__main__":
    main()