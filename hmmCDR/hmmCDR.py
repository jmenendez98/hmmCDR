import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os

from hmmlearn import hmm

from hmmCDR_parser import hmmCDR_parser
from hmmCDR_priors import hmmCDR_priors


class hmmCDR:
    '''
    CLASS DOCSTRING
    '''
    def __init__(self, bedgraphMethyl, output_path,
                 output_label, n_iter,
                 priors, emission_matrix, transition_matrix,
                 save_intermediates, use_percentiles, w, x, y, z):
        '''
        INIT DOCSTRING
        '''
        self.bedgraphMethyl = bedgraphMethyl
        self.output_label = output_label
        self.n_iter = n_iter
        self.output_path = output_path
        self.save_intermediates = save_intermediates
            
        self.use_percentiles = use_percentiles
        if self.use_percentiles is True:
            self.w = self.calculate_percentiles(self.bedgraphMethyl, w)
            self.x = self.calculate_percentiles(self.bedgraphMethyl, x)
            self.y = self.calculate_percentiles(self.bedgraphMethyl, y)
            self.z = self.calculate_percentiles(self.bedgraphMethyl, z)
        else:
            self.w = w
            self.x = x
            self.y = y
            self.z = z
        
        if priors is not None:
            # Check the type of priors and handle accordingly
            if isinstance(priors, pd.DataFrame):
                self.priors_df = priors
            elif isinstance(priors, str):
                self.priors_df = pd.read_csv(priors, sep='\t', header=None)
            else:
                raise ValueError("priors must be a pandas DataFrame or a string path to a file.")
            
            labelled_bedgraphMethyl = self.assign_emmisions(self.assign_priors(self.bedgraphMethyl))

            self.emission_matrix = self.calculate_emission_matrix(labelled_bedgraphMethyl)
            self.transition_matrix = self.calculate_transition_matrix(labelled_bedgraphMethyl)
        elif emission_matrix is not None and transition_matrix is not None:
            labelled_bedgraphMethyl = self.assign_emmisions(self.bedgraphMethyl)

            self.emission_matrix = emission_matrix
            self.transition_matrix = transition_matrix
            self.priors = None
        else:
            raise ValueError("Invalid input: Provide either a priors file or both emission and transition matrix files.")
            
        self.labelled_bedgraphMethyl = self.runHMM(labelled_bedgraphMethyl, 
                                                   transition_matrix=self.transition_matrix, 
                                                   emission_matrix=self.emission_matrix)
        self.hmmCDR_df = self.create_hmmCDR_df(self.labelled_bedgraphMethyl)
        self.hmmCDR_df.to_csv(self.output_path, sep="\t", header=None)

        
    def calculate_percentiles(self, bedgraphMethyl, percentile):
        '''
        DOCSTRING
        '''
        bedgraphMethyl['name'].replace('.', np.nan, inplace=True)
        methylation_scores = bedgraphMethyl['name'].dropna().astype(float)
        value = np.percentile(methylation_scores, q=percentile)
        return value
    
    def assign_priors(self, bedgraphMethyl):
        '''
        DOCSTRING
        '''
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bedgraphMethyl) # Convert DataFrames to BedTool objects
        priorCDR_df = self.priors_df[self.priors_df[3] == f'{self.output_label}'].drop(columns=[3])
        priorCDR_bedtool = pybedtools.BedTool.from_dataframe(priorCDR_df)
        priorTransition_df = self.priors_df[self.priors_df[3] == f'{self.output_label}_transition'].drop(columns=[3])
        priorTransition_bedtool = pybedtools.BedTool.from_dataframe(priorTransition_df)
        bedgraphMethyl['state'] = 'A' # Create a column for CDR state, default to 'A' (Not in CDR)
        intersected_cdr = bedMethyl_bedtool.intersect(priorCDR_bedtool, wa=True, wb=True) # Label CDR regions as 'C'
        cdr_df = intersected_cdr.to_dataframe()
        bedgraphMethyl.loc[bedgraphMethyl['start'].isin(cdr_df['start']), 'state'] = 'C' # Update the state based on matching 'start' values
        intersected_transition = bedMethyl_bedtool.intersect(priorTransition_bedtool, wa=True, wb=True) # Label CDR transition regions as 'B'
        transition_df = intersected_transition.to_dataframe()
        bedgraphMethyl.loc[bedgraphMethyl['start'].isin(transition_df['start']), 'state'] = 'B' # Update the state based on matching 'start' values
        return bedgraphMethyl

    def assign_emmisions(self, bedgraphMethyl):
        '''
        DOCSTRING
        '''
        def emissions_helper(value):
            if value > self.z:
                return 3
            elif value > self.y:
                return 2
            elif value > self.x:
                return 1
            elif value >= self.w:
                return 0
        bedgraphMethyl['emission'] = bedgraphMethyl['name'].apply(emissions_helper)
        return bedgraphMethyl

    def calculate_transition_matrix(self, labeled_bedMethyl_df):
        '''
        DOCSTRING
        '''
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
        '''
        DOCSTRING
        '''
        emission_counts = labeled_bedMethyl.groupby(['state', 'emission']).size().unstack(fill_value=0) # Get the counts of each emission in each state
        emission_matrix = emission_counts.div(emission_counts.sum(axis=1), axis=0) # Normalize counts to probabilities to get the emission matrix
        return emission_matrix

    def runHMM(self, emission_labelled_bedgraphMethyl, transition_matrix, emission_matrix):
        '''
        Create an HMM model with specified emission and transition matrices. Then Fit the HMM model on the emission data.

        Parameters:
        emission_matrix (numpy.ndarray): The emission probabilities matrix.
        transition_matrix (numpy.ndarray): The transition probabilities matrix.
        n_states (int): Number of states in the HMM.

        Returns:

        numpy.ndarray: The predicted states.
        '''
        # Create an HMM instance with Multinomial emissions
        model = hmm.CategoricalHMM(n_components=3, n_features=4, n_iter=self.n_iter, init_params="")
        # Set the start probs
        model.startprob_ = np.array([1.0, 0.0, 0.0])
        # Set the transition matrix
        model.transmat_ = transition_matrix
        #print('TransMat:', model.transmat_.shape, '\n', model.transmat_)
        # Set the emission matrix
        model.emissionprob_ = emission_matrix
        #print('EmisMat:', model.emissionprob_.shape, '\n', model.emissionprob_)
        emission_data = emission_labelled_bedgraphMethyl['emission'].to_numpy().reshape(-1, 1)
        # Predict the hidden states
        logprob, predicted_states = model.decode(emission_data, algorithm="viterbi")
        emission_labelled_bedgraphMethyl['predicted_state'] = predicted_states
        return emission_labelled_bedgraphMethyl

    def create_hmmCDR_df(self, df):
        '''
        DOCSTRING
        '''
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
        combined_df = pd.concat([merged_cdr_df, merged_transition_df], ignore_index=True) # Combine the two DataFrames
        combined_df = combined_df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
        def fix_transitions(df, distance=1000):
            df = df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
            for i in range(len(df)):
                if df.iloc[i]['name'] == 'CDR_transition':
                    prev_row = df.iloc[i - 1] if i > 0 else None
                    next_row = df.iloc[i + 1] if i < len(df) - 1 else None
                    if prev_row is not None and prev_row['name'] == 'CDR':
                        if df.iloc[i]['start'] <= prev_row['end'] + distance: # Check distance from transition start to CDR end
                            df.at[i, 'start'] = prev_row['end'] + 1 # Adjust start and end to be adjacent
                    if next_row is not None and next_row['name'] == 'CDR':
                        if next_row['start'] <= df.iloc[i]['end'] + distance: # Check distance from transition end to CDR start
                            df.at[i, 'end'] = next_row['start'] - 1 # Adjust start and end to be adjacent
            df = df[df['start'] <= df['end']] # Remove any transitions that end before they start (possible due to adjustment)
            return df
        adjusted_df = fix_transitions(combined_df, 1000)
        return adjusted_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process input files with optional parameters.')

    # Required arguments
    parser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    parser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    parser.add_argument('output_prefix', type=str, help='Output Prefix for the output files')

    # Parser Flags
    parser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    parser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Satellite type/name to filter CenSat bed file. (default: "H1L")')
    parser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')

    # Priors Flags
    parser.add_argument('--window_size', type=int, default=1020, help='Window size to calculate prior regions. (default: 1020)')
    parser.add_argument('--priorCDR_percent', type=int, default=5, help='Percentile for finding priorCDR regions. (default: 5)')
    parser.add_argument('--priorTransition_percent', type=int, default=10, help='Percentile for finding priorTransition regions. (default: 10)')
    parser.add_argument('--minCDR_size', type=int, default=3000, help='Minimum size for CDR regions. (default: 3000)')
    parser.add_argument('--enrichment', action='store_true', help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')

    # Create mutually exclusive group for file input options
    inputs_group = parser.add_mutually_exclusive_group(required=False)
    # Option for providing a single BED file
    inputs_group.add_argument('--cdr_priors', default=None, type=str, help='Path to the priorCDR bedfile')
    # Option for providing two matrix TSV files
    inputs_group.add_argument('--emission_matrix', default=None, type=str, help='Path to the emission matrix TSV file')
    inputs_group.add_argument('--transition_matrix', default=None, type=str, help='Path to the transition matrix TSV file')

    # HMM Flags
    parser.add_argument('--use_percentiles', action='store_true', default=True, help='Use values for flags w,x,y,z as percentile cutoffs for each category. (default: False)')
    parser.add_argument('--n_iter', type=int, default=1, help='Maximum number of iteration allowed for the HMM. (default: 1)')
    parser.add_argument('-w', type=int, default=0, help='Theshold for methylation to be classified as very low (default: 0)')
    parser.add_argument('-x', type=int, default=25, help='Theshold for methylation to be classified as low (default: 25)')
    parser.add_argument('-y', type=int, default=50, help='Theshold for methylation to be classified as medium (default: 50)')
    parser.add_argument('-z', type=int, default=75, help='Theshold for methylation to be classified as high (default: 75)')

    # Shared Flags
    parser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    parser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "CDR")')

    args = parser.parse_args()

    CDRparser = hmmCDR_parser(bedMethyl_path=args.bedMethyl_path,
                              cenSat_path=args.cenSat_path,
                              mod_code=args.mod_code,
                              sat_type=args.sat_type,
                              bedgraph=args.bedgraph
    )

    if not args.cdr_priors:
        CDRpriors = hmmCDR_priors(bedgraphMethyl=CDRparser.subset_bedgraphMethyl,
                                output_path=f'{args.output_prefix}_priorCDR.bed',
                                window_size=args.window_size, 
                                minCDR_size=args.minCDR_size, 
                                priorCDR_percent=args.priorCDR_percent, 
                                priorTransition_percent=args.priorTransition_percent, 
                                enrichment=args.enrichment, 
                                save_intermediates=args.save_intermediates, 
                                output_label=args.output_label
        )

        CDRhmm = hmmCDR(bedgraphMethyl=CDRparser.subset_bedgraphMethyl,
                        priors=CDRpriors.hmmCDRpriors,
                        output_path=f'{args.output_prefix}_hmmCDR.bed',
                        output_label=args.output_label,
                        n_iter=args.n_iter,
                        emission_matrix=args.emission_matrix,
                        transition_matrix=args.transition_matrix,
                        use_percentiles=args.use_percentiles,
                        save_intermediates=args.save_intermediates,
                        w=args.w, x=args.x, y=args.y, z=args.z
        )
    else:
        CDRhmm = hmmCDR(bedgraphMethyl=CDRparser.subset_bedgraphMethyl,
                        priors=args.cdr_priors,
                        output_path=args.output_path,
                        output_label=args.output_label,
                        n_iter=args.n_iter,
                        emission_matrix=args.emission_matrix,
                        transition_matrix=args.transition_matrix,
                        use_percentiles=args.use_percentiles,
                        w=args.w, x=args.x, y=args.y, z=args.z
        )