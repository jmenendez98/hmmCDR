import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os
import concurrent.futures

from hmmlearn import hmm

from hmmCDR.parser import hmmCDR_parser
from hmmCDR.find_priors import hmmCDR_prior_finder


class hmmCDR:
    def __init__(self, output_label, n_iter, min_size, merge_distance, 
                 main_color, transition_color, remove_transitions, raw_thresholds,
                 emission_matrix=None, transition_matrix=None, w=0, x=25, y=50, z=75):

        self.output_label = output_label
        self.n_iter = n_iter
        self.min_size = min_size
        self.merge_distance = merge_distance
        self.raw_thresholds = raw_thresholds

        self.remove_transitions = remove_transitions
        self.main_color = main_color
        self.transition_color = transition_color 

        self.w, self.x, self.y, self.z = w, x, y, z

        self.emission_matrix = emission_matrix
        self.transition_matrix = transition_matrix

    
    def assign_priors(self, bed4Methyl, priors):
        bed4Methyl['prior'] = 'A'

        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bed4Methyl)
        prior_bedtool = pybedtools.BedTool.from_dataframe(priors)

        intersected_df = bedMethyl_bedtool.intersect(prior_bedtool, wa=True, wb=True).to_dataframe()
        bed4Methyl.loc[bed4Methyl['start'].isin(intersected_df['start']), 'prior'] = 'C'

        return bed4Methyl

    def calculate_transition_matrix(self, labeled_bedMethyl_df):
        transitions = {f'{a}->{b}': 0 for a in 'ABC' for b in 'ABC'}
        
        prev_state = labeled_bedMethyl_df.iloc[0]['prior']
        for state in labeled_bedMethyl_df['prior'][1:]:
            transitions[f'{prev_state}->{state}'] += 1
            prev_state = state
        totals = {state: sum(transitions[f'{state}->{next_state}'] for next_state in 'ABC') for state in 'ABC'}
        
        transition_matrix = [[transitions[f'{a}->{b}'] / totals[a] if totals[a] else 0 for b in 'ABC'] for a in 'ABC']

        initial_A_C = transition_matrix[0][2]  # A->C value before modification
        transition_matrix[0][1] = initial_A_C / 3  # A->B becomes 1/3 of A->C
        transition_matrix[0][2] = 2 * initial_A_C / 3  # A->C becomes 2/3 of its original value

        initial_C_A = transition_matrix[2][0]  # C->A value before modification
        transition_matrix[2][0] = initial_C_A / 3  # C->A becomes 1/3 of its original value
        transition_matrix[2][1] = 2 * initial_C_A / 3  # C->B becomes 2/3 of the initial C->A

        transition_matrix[1][0] = initial_C_A / 3  # B->A becomes 1/3 of initial C->A
        transition_matrix[1][2] = 2 * initial_C_A / 3  # B->C becomes 2/3 of initial C->A
        transition_matrix[1][1] = transition_matrix[2][2]  # B->B becomes the same as C->C

        return transition_matrix
    
    def calculate_emission_thresholds(self, bed4Methyl):
        if self.raw_thresholds:
            return sorted([self.w, self.x, self.y, self.z])
        
        methylation_scores = pd.to_numeric(bed4Methyl['name'].replace('.', np.nan), errors='coerce').dropna()
        methylation_scores = [0] + methylation_scores[methylation_scores != 0].tolist()
        
        return sorted(np.percentile(methylation_scores, q=[self.w, self.x, self.y, self.z]))

    def assign_emissions(self, bed4Methyl, emission_thresholds):
        def emissions_helper(value):
            for i, threshold in enumerate(emission_thresholds):
                if value <= threshold:
                    return i
            return len(emission_thresholds)-1    
        
        bed4Methyl['emission'] = bed4Methyl['name'].apply(emissions_helper)
        return bed4Methyl

    def calculate_emission_matrix(self, labeled_bedMethyl):
        state_mapping = {'A': 0, 'B': 1, 'C': 2}

        emission_matrix = np.zeros((3, 4))
        emission_counts = labeled_bedMethyl.groupby(['prior', 'emission']).size().unstack(fill_value=0)

        for prior, i in state_mapping.items():
            if prior not in emission_counts.index:
                new_row = pd.Series([0] * len(emission_counts.columns), index=emission_counts.columns, name=prior)
                # Instead of using append, use pd.concat
                emission_counts = pd.concat([emission_counts, new_row.to_frame().T])

        for emission in range(4):
            if emission not in emission_counts.columns:
                emission_counts[emission] = 0

        # Fill emission matrix for A and C
        for prior, i in state_mapping.items():
            if prior in emission_counts.index and prior != 'B':  # Handle A and C first
                total = emission_counts.loc[prior].sum()

                for j in range(4):
                    if j in emission_counts.columns and total > 0:
                        emission_matrix[i, j] = emission_counts.loc[prior, j] / total
                    else:
                        emission_matrix[i, j] = 0

        # Compute B row as the average of the A and C rows
        emission_matrix[1, :] = (emission_matrix[0, :] + emission_matrix[2, :]) / 2
        
        # Normalize the B row so that it sums to 1
        b_row_sum = emission_matrix[1, :].sum()
        if b_row_sum > 0:
            emission_matrix[1, :] = emission_matrix[1, :] / b_row_sum

        return emission_matrix


    def runHMM(self, emission_labelled_bed4Methyl, transition_matrix, emission_matrix):
        model = hmm.CategoricalHMM(n_components=3, n_iter=self.n_iter, init_params="")
        model.startprob_ = np.array([1.0, 0.0, 0.0])
        model.transmat_ = transition_matrix
        model.emissionprob_ = emission_matrix

        emission_data = emission_labelled_bed4Methyl['emission'].values.reshape(-1, 1)
        _, predicted_states = model.decode(emission_data, algorithm="viterbi")

        emission_labelled_bed4Methyl['predicted_state'] = predicted_states
        return emission_labelled_bed4Methyl

    def create_hmmCDR_df(self, df):
        '''
        DOCSTRING
        '''
        def merge_and_label_regions(df, state_value, label, color):
            state_df = df[df['predicted_state'] == state_value][['chrom', 'start', 'end']]
            if not state_df.empty:
                state_bedtool = pybedtools.BedTool.from_dataframe(state_df)
                merged_df = state_bedtool.merge(d=self.merge_distance).to_dataframe(names=['chrom', 'start', 'end'])
                merged_df['name'] = label
                merged_df = merged_df[(merged_df['end'] - merged_df['start']) >= self.min_size]
                merged_df['score'] = 0
                merged_df['strand'] = '.'
                merged_df['thickStart'] = merged_df['start']
                merged_df['thickEnd'] = merged_df['end']
                merged_df['color'] = color
                return merged_df
            else:
                return pd.DataFrame()
        
        # Merge and label CDR and transition regions
        merged_cdr_df = merge_and_label_regions(df, 2, f'{self.output_label}', self.main_color)
        merged_transition_df = merge_and_label_regions(df, 1, f'{self.output_label}_transition', self.transition_color)

        if not self.remove_transitions:
            return merged_cdr_df

        combined_df = pd.concat([merged_cdr_df, merged_transition_df]).sort_values(by=['chrom', 'start']).reset_index(drop=True)

        def fix_transitions(df):
            for i in range(len(df)):
                if df.iloc[i]['name'] == f'{self.output_label}_transition':
                    prev_row = df.iloc[i - 1] if i > 0 else None
                    next_row = df.iloc[i + 1] if i < len(df) - 1 else None
                    if prev_row is not None and prev_row['name'] == self.output_label:
                        if df.iloc[i]['start'] <= prev_row['end'] + self.merge_distance: 
                            df.at[i, 'start'] = prev_row['end'] + 1 
                    if next_row is not None and next_row['name'] == self.output_label:
                        if next_row['start'] <= df.iloc[i]['end'] + self.merge_distance: 
                            df.at[i, 'end'] = next_row['start'] - 1 
            return df[df['start'] <= df['end']]

        return fix_transitions(combined_df)
    
    def hmm_single_chromosome(self, chrom, bed4Methyl_chrom, priors_chrom):
        labelled_bed4Methyl_chrom = self.assign_emissions(
            self.assign_priors(bed4Methyl_chrom, priors_chrom), 
            self.calculate_emission_thresholds(bed4Methyl_chrom)
        )

        emission_matrix = self.emission_matrix or self.calculate_emission_matrix(labelled_bed4Methyl_chrom)
        transition_matrix = self.transition_matrix or self.calculate_transition_matrix(labelled_bed4Methyl_chrom)

        hmmlabelled_bed4Methyl = self.runHMM(labelled_bed4Methyl_chrom, transition_matrix, emission_matrix)
        single_chrom_hmmCDR_results = self.create_hmmCDR_df(hmmlabelled_bed4Methyl)

        return chrom, single_chrom_hmmCDR_results, hmmlabelled_bed4Methyl

    def hmm_all_chromosomes(self, bed4Methyl_chrom_dict, priors_chrom_dict):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(self.hmm_single_chromosome, chrom, bed4Methyl_chrom_dict[chrom], priors_chrom_dict[chrom]): chrom 
                for chrom in bed4Methyl_chrom_dict
            }

            results = {chrom: future.result() for future, chrom in futures.items()}
            hmmCDRresults_chrom_dict = {chrom: result[1] for chrom, result in results.items()}
            hmmCDR_labelled_bed4Methyl_chrom_dict = {chrom: result[2] for chrom, result in results.items()}

        self.chromosomes = list(bed4Methyl_chrom_dict.keys())
        self.hmmCDRpriors_chrom_dict = hmmCDRresults_chrom_dict

        return hmmCDRresults_chrom_dict, hmmCDR_labelled_bed4Methyl_chrom_dict

def main():
    argparser= argparse.ArgumentParser(description='Process input files with optional parameters.')

    argparser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    argparser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    argparser.add_argument('output_path', type=str, help='Output Path for the output files')

    # hmmCDR Parser Flags
    argparser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    argparser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")')
    argparser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--min_valid_cov', type=int, default=10, help='Minimum Valid Coverage to consider a methylation site. (default: 10)')

    # hmmCDR Priors Flags
    argparser.add_argument('--window_size', type=int, default=510, help='Window size to calculate prior regions. (default: 510)')
    argparser.add_argument('--prior_percentile', type=float, default=10, help='Percentile for finding priorCDR regions. (default: 10)')

    # HMM Flags
    argparser.add_argument('--raw_thresholds', action='store_true', default=False, help='Use values for flags w,x,y,z as raw threshold cutoffs for each emission category. (default: True)')
    argparser.add_argument('--n_iter', type=int, default=1, help='Maximum number of iteration allowed for the HMM. (default: 1)')
    argparser.add_argument('--remove_transitions', action='store_true', default=False, help='Do not report transitions in final hmmCDR output file. (default: False)')
    argparser.add_argument('-w', type=int, default=0, help='Threshold of non-zero methylation percentile to be classified as very low (default: 0)')
    argparser.add_argument('-x', type=int, default=25, help='Threshold of non-zero methylation percentile to be classified as low (default: 25)')
    argparser.add_argument('-y', type=int, default=50, help='Threshold of non-zero methylation percentile to be classified as medium (default: 50)')
    argparser.add_argument('-z', type=int, default=75, help='Threshold of non-zero methylation percentile to be classified as high (default: 75)')

    # Shared Flags
    argparser.add_argument('--merge_distance', type=int, default=1021, help='Distance to merge adjacently labelled regions. (default: 1021)')
    argparser.add_argument('--min_size', type=int, default=3000, help='Minimum size for regions. (default: 3000)')
    argparser.add_argument('--enrichment', action='store_true', default=False, help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')

    argparser.add_argument('--main_color', type=str, default='50,50,255', help='Color to dictate main regions. (default: 50,50,255)')
    argparser.add_argument('--transition_color', type=str, default='100,150,200', help='Color to dictate transition regions. (default: 100,150,200)')
    argparser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    argparser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "CDR")')

    args = argparser.parse_args()
    output_prefix = os.path.splitext(args.output_path)[0]
    sat_types = [st.strip() for st in args.sat_type.split(',')]

    CDRparser = hmmCDR_parser(
        mod_code=args.mod_code,
        sat_type=sat_types,
        bedgraph=args.bedgraph,
        min_valid_cov=args.min_valid_cov
    )

    bed4Methyl_chrom_dict, cenSat_chrom_dict = CDRparser.process_files(
        bedMethyl_path=args.bedMethyl_path, 
        cenSat_path=args.cenSat_path
    )

    if args.save_intermediates:
        concatenated_bed4Methyl = pd.concat(bed4Methyl_chrom_dict.values(), axis=0)
        bed4Methyl_output_path = f'{output_prefix}_bed4.bedgraph'
        concatenated_bed4Methyl.to_csv(bed4Methyl_output_path, sep='\t', index=False, header=False)
        print(f'Wrote methylation bedgraph to: {bed4Methyl_output_path}')
    
    CDRpriors = hmmCDR_prior_finder(
        window_size=args.window_size, 
        min_size=args.min_size, 
        prior_percentile=args.prior_percentile, 
        merge_distance=args.merge_distance, 
        enrichment=args.enrichment, 
        output_label=args.output_label
    )

    hmmCDRpriors_chrom_dict = CDRpriors.priors_all_chromosomes(bed4Methyl_chrom_dict=bed4Methyl_chrom_dict)
    
    if args.save_intermediates:
        concatenated_priors = pd.concat(hmmCDRpriors_chrom_dict.values(), axis=0)
        intermediate_output_path = f'{output_prefix}_priors.bed'
        concatenated_priors.to_csv(intermediate_output_path, sep='\t', index=False, header=False)
        print(f'Wrote priors to: {intermediate_output_path}')

    CDRhmm = hmmCDR(
        output_label=args.output_label,
        n_iter=args.n_iter,
        min_size=args.min_size,
        merge_distance=args.merge_distance, 
        remove_transitions=args.remove_transitions,
        main_color=args.main_color,
        transition_color=args.transition_color,
        raw_thresholds=args.raw_thresholds,
        w=args.w, x=args.x, y=args.y, z=args.z
    )

    hmmCDRresults_chrom_dict, hmm_labelled_bed4Methyl_chrom_dict = CDRhmm.hmm_all_chromosomes(
        bed4Methyl_chrom_dict=bed4Methyl_chrom_dict,
        priors_chrom_dict=hmmCDRpriors_chrom_dict
    )

    concatenated_hmmCDRs = pd.concat(hmmCDRresults_chrom_dict.values(), axis=0)
    concatenated_hmmCDRs.to_csv(args.output_path, sep='\t', index=False, header=False)
    if args.save_intermediates: 
        print(f"Final output saved to: {args.output_path}")


if __name__ == "__main__":
    main()
    