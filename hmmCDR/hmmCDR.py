import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os
import concurrent.futures

from hmmlearn import hmm

from hmmCDR.hmmCDRparse import hmmCDRparse
from hmmCDR.hmmCDRprior import hmmCDRprior


class hmmCDR:
    '''
    CLASS DOCSTRING
    '''
    def __init__(self, output_label, n_iter,
                 emission_matrix=None, transition_matrix=None,
                 use_percentiles=False,
                 w=0, x=25, y=50, z=75):
        '''
        INIT DOCSTRING
        '''
        # all hmmCDR class parameters are optional
        self.output_label = output_label
        self.n_iter = n_iter
        self.use_percentiles = use_percentiles

        self.w, self.x, self.y, self.z = w, x, y, z

        self.emission_matrix = emission_matrix
        self.transition_matrix = transition_matrix
        '''
        if priors is not None:
            # Check the type of priors and handle accordingly
            if isinstance(priors, pd.DataFrame):
                priors = priors
            elif isinstance(priors, str):
                priors = pd.read_csv(priors, sep='\t', header=None)
            else:
                raise ValueError("priors must be a pandas DataFrame or a string path to a file.")
            
            labelled_bed4Methyl = self.assign_emmisions(self.assign_priors(self.bed4Methyl))

            self.emission_matrix = self.calculate_emission_matrix(labelled_bed4Methyl)
            self.transition_matrix = self.calculate_transition_matrix(labelled_bed4Methyl)
        elif emission_matrix is not None and transition_matrix is not None:
            labelled_bed4Methyl = self.assign_emmisions(self.bed4Methyl)

            self.emission_matrix = emission_matrix
            self.transition_matrix = transition_matrix
            self.priors = None
        else:
            raise ValueError("Invalid input: Provide either a priors file or both emission and transition matrix files.")
            
        self.labelled_bed4Methyl = self.runHMM(labelled_bed4Methyl, 
                                                   transition_matrix=self.transition_matrix, 
                                                   emission_matrix=self.emission_matrix)
        self.hmmCDR_df = self.create_hmmCDR_df(self.labelled_bed4Methyl)
        self.hmmCDR_df.to_csv(self.output_path, sep="\t", header=None)        
        '''
    
    def assign_priors(self, bed4Methyl, priors):
        '''
        DOCSTRING
        '''
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bed4Methyl) # Convert DataFrames to BedTool objects
        priorCDR_df = priors[priors[3] == f'{self.output_label}'].drop(columns=[3])
        priorCDR_bedtool = pybedtools.BedTool.from_dataframe(priorCDR_df)
        priorTransition_df = priors[priors[3] == f'{self.output_label}_transition'].drop(columns=[3])
        priorTransition_bedtool = pybedtools.BedTool.from_dataframe(priorTransition_df)
        bed4Methyl['state'] = 'A' # Create a column for CDR state, default to 'A' (Not in CDR)
        intersected_cdr = bedMethyl_bedtool.intersect(priorCDR_bedtool, wa=True, wb=True) # Label CDR regions as 'C'
        cdr_df = intersected_cdr.to_dataframe()
        bed4Methyl.loc[bed4Methyl['start'].isin(cdr_df['start']), 'state'] = 'C' # Update the state based on matching 'start' values
        intersected_transition = bedMethyl_bedtool.intersect(priorTransition_bedtool, wa=True, wb=True) # Label CDR transition regions as 'B'
        transition_df = intersected_transition.to_dataframe()
        bed4Methyl.loc[bed4Methyl['start'].isin(transition_df['start']), 'state'] = 'B' # Update the state based on matching 'start' values
        return bed4Methyl

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
    
    def calculate_emission_thresholds(self, bed4Methyl):
        bed4Methyl['name'].replace('.', np.nan, inplace=True)
        methylation_scores = bed4Methyl['name'].dropna().astype(float)
        if self.use_percentiles:
            return [np.percentile(methylation_scores, q=percentile) for percentile in [self.w, self.x, self.y, self.z]]
        else:
            return [self.w, self.x, self.y, self.z]

    def assign_emmisions(self, bed4Methyl, emission_thesholds):
        '''
        DOCSTRING
        '''
        def emissions_helper(value):
            if value > emission_thesholds[3]:
                return 3
            elif value > emission_thesholds[2]:
                return 2
            elif value > emission_thesholds[1]:
                return 1
            elif value >= emission_thesholds[0]:
                return 0
        bed4Methyl['emission'] = bed4Methyl['name'].apply(emissions_helper)
        return bed4Methyl

    def calculate_emission_matrix(self, labeled_bedMethyl):
        '''
        DOCSTRING
        '''
        emission_counts = labeled_bedMethyl.groupby(['state', 'emission']).size().unstack(fill_value=0) # Get the counts of each emission in each state
        emission_matrix = emission_counts.div(emission_counts.sum(axis=1), axis=0) # Normalize counts to probabilities to get the emission matrix
        return emission_matrix

    def runHMM(self, emission_labelled_bed4Methyl, transition_matrix, emission_matrix):
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
        emission_data = emission_labelled_bed4Methyl['emission'].to_numpy().reshape(-1, 1)
        # Predict the hidden states
        logprob, predicted_states = model.decode(emission_data, algorithm="viterbi")
        emission_labelled_bed4Methyl['predicted_state'] = predicted_states
        return emission_labelled_bed4Methyl

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
    
    def hmm_single_chromosome(self, chrom, bed4Methyl_chrom, priors_chrom):
        if self.emission_matrix is None and self.transition_matrix is None:
            labelled_bed4Methyl_chrom = self.assign_emmisions(self.assign_priors(bed4Methyl_chrom, priors_chrom), 
                                                              self.calculate_emission_thresholds(bed4Methyl_chrom))
            emission_matrix = self.calculate_emission_matrix(labelled_bed4Methyl_chrom)
            transition_matrix = self.calculate_transition_matrix(labelled_bed4Methyl_chrom)
        else:
            labelled_bed4Methyl_chrom = self.assign_emmisions(bed4Methyl_chrom, 
                                                              self.calculate_emission_thresholds(bed4Methyl_chrom))
            emission_matrix = self.emission_matrix
            transition_matrix = self.transition_matrix
        hmmlabelled_bed4Methyl = self.runHMM(labelled_bed4Methyl_chrom,
                                             transition_matrix=transition_matrix,
                                             emission_matrix=emission_matrix)
        single_chrom_hmmCDR_results = self.create_hmmCDR_df(hmmlabelled_bed4Methyl)
        return chrom, single_chrom_hmmCDR_results, hmmlabelled_bed4Methyl

    def hmm_all_chromosomes(self, bed4Methyl_chrom_dict, priors_chrom_dict):
        '''
        Processes all chromosomes in parallel using concurrent futures.

        Parameters:
        -----------
        bed4Methyl_chrom_dict : dict
            A dictionary with chromosome names as keys and DataFrames (from bedMethyl) as values.
        
        cenSat : pd.DataFrame
            The DataFrame containing cenSat annotations.

        Returns:
        --------
        dict
            Two dictionaries:
            - final_bed4Methyl_chrom_dict: A dictionary with chromosome names as keys and processed DataFrames as values.
            - cenSat_chrom_dict: A dictionary with chromosome names as keys and filtered cenSat DataFrames as values.
        '''
        hmmCDRresults_chrom_dict = {}
        hmmCDR_labelled_bed4Methyl_chrom_dict = {}
        chromosomes = bed4Methyl_chrom_dict.keys()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.hmm_single_chromosome, chrom,
                    bed4Methyl_chrom_dict[chrom],
                    priors_chrom_dict[chrom]
                ): chrom for chrom in chromosomes
            }

            for future in concurrent.futures.as_completed(futures):
                chrom, hmmCDR_chrom_result, hmmCDR_labelled_bed4Methyl = future.result()

                hmmCDRresults_chrom_dict[chrom] = hmmCDR_chrom_result
                hmmCDR_labelled_bed4Methyl_chrom_dict[chrom] = hmmCDR_labelled_bed4Methyl

        self.chromosomes = chromosomes
        self.hmmCDRpriors_chrom_dict = hmmCDRresults_chrom_dict

        return hmmCDRresults_chrom_dict, hmmCDR_labelled_bed4Methyl_chrom_dict

def main():
    argparser= argparse.ArgumentParser(description='Process input files with optional parameters.')

    # If --matrix is passed in the first two inputs are paths to matrix files.
    argparser.add_argument('--matrix', action='store_true', help='Use matrices instead of BED files. (default: False)')
    required_args = argparser.add_argument_group('required arguments')
    matrix_arg, unknown = argparser.parse_known_args()

    if matrix_arg.matrix:
        required_args.add_argument('transition_matrix', type=str, help='Transition Matrix')
        required_args.add_argument('emission_matrix', type=str, help='Emission Matrix')
    else:
        required_args.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
        required_args.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    argparser.add_argument('output_path', type=str, help='Output Path for the output files')

    # hmmCDR Parser Flags
    argparser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--rolling_window', type=int, default=0, help='Flag indicating whether or not to use a rolling average and the rolling avg window size. If set to 0 no rolling averages are used. (defualt: 0)')
    argparser.add_argument('--min_valid_cov', type=int, default=10, help='Minimum Valid Coverage to consider a methylation site. (default: 10)')
    argparser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    argparser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Satellite type/name to filter CenSat bed file. (default: "H1L")')

    # hmmCDR Priors Flags
    argparser.add_argument('--window_size', type=int, default=1020, help='Window size to calculate prior regions. (default: 1020)')
    argparser.add_argument('--priorCDR_percent', type=int, default=5, help='Percentile for finding priorCDR regions. (default: 5)')
    argparser.add_argument('--priorTransition_percent', type=int, default=10, help='Percentile for finding priorTransition regions. (default: 10)')
    argparser.add_argument('--minCDR_size', type=int, default=3000, help='Minimum size for CDR regions. (default: 3000)')
    argparser.add_argument('--enrichment', action='store_true', help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    argparser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")

    # HMM Flags
    argparser.add_argument('--use_percentiles', action='store_true', default=False, help='Use values for flags w,x,y,z as percentile cutoffs for each category. (default: False)')
    argparser.add_argument('--n_iter', type=int, default=1, help='Maximum number of iteration allowed for the HMM. (default: 1)')
    argparser.add_argument('-w', type=int, default=0, help='Theshold for methylation to be classified as very low (default: 0)')
    argparser.add_argument('-x', type=int, default=25, help='Theshold for methylation to be classified as low (default: 25)')
    argparser.add_argument('-y', type=int, default=50, help='Theshold for methylation to be classified as medium (default: 50)')
    argparser.add_argument('-z', type=int, default=75, help='Theshold for methylation to be classified as high (default: 75)')

    # Shared Flags
    argparser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    argparser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "CDR")')

    args = argparser.parse_args()
    output_prefix = os.path.splitext(args.output_path)[0]

    # Extract required arguments as variables
    transition_matrix = getattr(args, 'transition_matrix', None)
    emission_matrix = getattr(args, 'emission_matrix', None)
    bedMethyl_path = getattr(args, 'bedMethyl_path', None)
    cenSat_path = getattr(args, 'cenSat_path', None)

    if not args.matrix:
        CDRparser = hmmCDRparse(
            bedMethyl_path=bedMethyl_path,
            cenSat_path=cenSat_path,
            mod_code=args.mod_code,
            sat_type=args.sat_type,
            bedgraph=args.bedgraph,
            min_valid_cov=args.min_valid_cov,
            rolling_window=args.rolling_window
        )

        cenSat = CDRparser.read_cenSat(path=CDRparser.cenSat_path)
        bedMethyl = CDRparser.read_bedMethyl(path=CDRparser.bedMethyl_path)
        bed4Methyl_chrom_dict, cenSat_chrom_dict = CDRparser.parse_all_chromosomes(bedMethyl=bedMethyl, cenSat=cenSat)

        if args.save_intermediates:
            concatenated_bed4Methyls = pd.concat(bed4Methyl_chrom_dict.values(), axis=0)
            concatenated_bed4Methyls.to_csv(f'{output_prefix}_intersected_bed4Methyl.bedgraph', 
                                            sep='\t', index=False, header=False)
            concatenated_regions = pd.concat(cenSat_chrom_dict.values(), axis=0)
            concatenated_regions.to_csv(f'{output_prefix}_selected_regions.bed', 
                                            sep='\t', index=False, header=False)
            print(f'Wrote Intermediates: {output_prefix}_intersected_bed4Methyl.bedgraph and {output_prefix}_selected_regions.bed.')

        CDRpriors = hmmCDRprior(
            window_size=args.window_size, 
            minCDR_size=args.minCDR_size, 
            priorCDR_percent=args.priorCDR_percent, 
            priorTransition_percent=args.priorTransition_percent, 
            enrichment=args.enrichment, 
            output_label=args.output_label
        )

        hmmCDRpriors_chrom_dict = CDRpriors.priors_all_chromosomes(bed4Methyl_chrom_dict=bed4Methyl_chrom_dict)
        
        if args.save_intermediates:
            concatenated_priors = pd.concat(hmmCDRpriors_chrom_dict.values(), axis=0)
            concatenated_priors.to_csv(f'{output_prefix}_hmmCDRpriors.bed', 
                                            sep='\t', index=False, header=False)
            print(f'Wrote Intermediate: {output_prefix}_hmmCDRpriors.bed.')

    CDRhmm = hmmCDR(
        output_label=args.output_label,
        n_iter=args.n_iter,
        emission_matrix=emission_matrix,
        transition_matrix=transition_matrix,
        use_percentiles=args.use_percentiles,
        w=args.w, x=args.x, y=args.y, z=args.z
    )

    hmmCDRresults_chrom_dict, hmm_labelled_bed4Methyl_chrom_dict = CDRhmm.hmm_all_chromosomes(bed4Methyl_chrom_dict=bed4Methyl_chrom_dict, 
                                                                                              priors_chrom_dict=hmmCDRpriors_chrom_dict)

    if args.save_intermediates:
        concatenated_hmmCDR_labelled_bed4Methyl = pd.concat(hmm_labelled_bed4Methyl_chrom_dict.values(), axis=0)
        concatenated_hmmCDR_labelled_bed4Methyl.to_csv(f'{output_prefix}_hmmCDR_labelled_bed4Methyl.bed', 
                                                       sep='\t', index=False, header=False)
        print(f'Wrote Intermediate: {output_prefix}_hmmCDR_labelled_bed4Methyl.bed.')

    # Combine all chromosomes and save the output
    concatenated_hmmCDRs = pd.concat(hmmCDRresults_chrom_dict.values(), axis=0)
    concatenated_hmmCDRs.to_csv(args.output_path, sep='\t', index=False, header=False)
    print(f"hmmCDRs saved to: {args.output_path}")


if __name__ == "__main__":
    main()
    