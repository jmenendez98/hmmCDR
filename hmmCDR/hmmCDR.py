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
    def __init__(self, raw_thresholds, n_iter, 
                 min_size, merge_distance, min_cdr_score, min_low_conf_score, 
                 main_color, low_conf_color, output_label,
                 w=0, x=25, y=50, z=75):

        self.raw_thresholds = raw_thresholds
        self.n_iter = n_iter
        
        self.min_size = min_size
        self.merge_distance = merge_distance
        self.min_cdr_score = min_cdr_score
        self.min_low_conf_score = min_low_conf_score

        self.main_color = main_color
        self.low_conf_color = low_conf_color 
        self.output_label = output_label

        self.w, self.x, self.y, self.z = w, x, y, z

    
    def assign_priors(self, bed4Methyl, priors):
        bed4Methyl['prior'] = 0

        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bed4Methyl)
        prior_bedtool = pybedtools.BedTool.from_dataframe(priors)

        intersected_df = bedMethyl_bedtool.intersect(prior_bedtool, wa=True, wb=True).to_dataframe()
        bed4Methyl.loc[bed4Methyl['start'].isin(intersected_df['start']), 'prior'] = 1

        return bed4Methyl

    def calculate_transition_matrix(self, labeled_bedMethyl_df):
        # Initialize the transitions dictionary for states '0' and '1'
        transitions = {'0->0': 0, '0->1': 0, '1->0': 0, '1->1': 0}
        
        # Track transitions between states
        prev_state = labeled_bedMethyl_df.iloc[0]['prior']
        for state in labeled_bedMethyl_df['prior'][1:]:
            transitions[f'{prev_state}->{state}'] += 1
            prev_state = state
        
        # Calculate totals for each state to normalize the transition counts
        totals = {state: transitions[f'{state}->0'] + transitions[f'{state}->1'] for state in ['0', '1']}
        
        # Create the transition matrix for states 0 and 1
        transition_matrix = [
            [
                transitions[f'{a}->{b}'] / totals[a] if totals[a] else 0
                for b in ['0', '1']
            ] for a in ['0', '1']
        ]
        
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
        state_mapping = {'0': 0, '1': 1}

        emission_matrix = np.zeros((2, 4))
        
        # Group the data by 'prior' (state) and 'emission', counting occurrences
        emission_counts = labeled_bedMethyl.groupby(['prior', 'emission']).size().unstack(fill_value=0)

        row_sums = emission_counts.sum(axis=1)

        emission_matrix = emission_counts.div(row_sums, axis=0)

        return emission_matrix

    def runHMM(self, emission_labelled_bed4Methyl, transition_matrix, emission_matrix):
        model = hmm.CategoricalHMM(n_components=2, n_iter=self.n_iter, init_params="")
        model.startprob_ = np.array([1.0, 0.0])
        model.transmat_ = transition_matrix
        model.emissionprob_ = emission_matrix

        emission_data = emission_labelled_bed4Methyl['emission'].values.reshape(-1, 1)
        _, predicted_states = model.decode(emission_data, algorithm="viterbi")
        log_likelihood, responsibilities = model.score_samples(emission_data)

        emission_labelled_bed4Methyl['CDR_score'] = [score[1] for score in responsibilities]

        return emission_labelled_bed4Methyl

    def create_subCDR_df(self, df):
        # isolate CpG positions and scores
        single_chrom_hmmCDR_scores = df[['chrom', 'start', 'end', 'CDR_score']]

        def merge_CpG_dataframe(dataframe, threshold):
            merged_rows = []
            score_list = []
            current_chrom, current_start, current_end = None, None, None
            dataframe = dataframe[dataframe['CDR_score'] > threshold]
            for _, row in dataframe.iterrows():
                chrom, start, end, score = row['chrom'], row['start'], row['end'], row['CDR_score']
                if current_chrom == chrom and current_end is not None and start - current_end <= self.merge_distance:
                    current_end = end
                    score_list.append(score)
                else:
                    if current_chrom is not None:
                        # Calculate the mean of scores and round to 5 decimals
                        avg_score = np.mean(score_list)
                        merged_rows.append([current_chrom, current_start, current_end, round(avg_score, 5)])
                    current_chrom, current_start, current_end = chrom, start, end
                    score_list = [score]
            if current_chrom is not None:
                avg_score = np.mean(score_list)
                merged_rows.append([current_chrom, current_start, current_end, round(avg_score, 5)])
            return pd.DataFrame(merged_rows, columns=['chrom', 'start', 'end', 'score'], index=None)

        # Create DataFrame for high confidence CDRs
        single_chrom_CDRs = merge_CpG_dataframe(single_chrom_hmmCDR_scores, self.min_cdr_score)
        single_chrom_CDRs['size'] = single_chrom_CDRs['end'] - single_chrom_CDRs['start']
        single_chrom_CDRs = single_chrom_CDRs[single_chrom_CDRs['size'] >= self.min_size]
        cdr_scores_column = single_chrom_CDRs['score']

        # Handle low confidence CDRs if thresholds are defined
        if self.min_low_conf_score > 0.0 and self.min_low_conf_score < self.min_cdr_score:
            single_chrom_low_conf_CDRs = merge_CpG_dataframe(single_chrom_hmmCDR_scores, self.min_low_conf_score)
            single_chrom_low_conf_CDRs['size'] = single_chrom_low_conf_CDRs['end'] - single_chrom_low_conf_CDRs['start']
            single_chrom_low_conf_CDRs = single_chrom_low_conf_CDRs[single_chrom_low_conf_CDRs['size'] >= self.min_size]

            # subtract high confidence from low confidence CDRs
            cdr_bedtool = pybedtools.BedTool.from_dataframe(single_chrom_CDRs)
            low_conf_cdr_bedtool = pybedtools.BedTool.from_dataframe(single_chrom_low_conf_CDRs)
            single_chrom_low_conf_CDRs = low_conf_cdr_bedtool.subtract(cdr_bedtool).to_dataframe(names=['chrom', 'start', 'end', 'score', 'size'])
            low_conf_scores_column = single_chrom_low_conf_CDRs['score']

            # Deal with if you removed all low confidence positions
            if not single_chrom_low_conf_CDRs.empty:
                single_chrom_low_conf_CDRs = single_chrom_low_conf_CDRs[['chrom', 'start', 'end']]
                single_chrom_low_conf_CDRs['name'] = f"low_conf_sub{self.output_label}"
                single_chrom_low_conf_CDRs['score'] = low_conf_scores_column * 100
                single_chrom_low_conf_CDRs['strand'] = '.'
        else:
            single_chrom_low_conf_CDRs = pd.DataFrame()

        # Filter and prepare the final CDR DataFrame
        single_chrom_CDRs = single_chrom_CDRs[['chrom', 'start', 'end']]
        single_chrom_CDRs['name'] = f"sub{self.output_label}"
        single_chrom_CDRs['score'] = cdr_scores_column * 100
        single_chrom_CDRs['strand'] = '.'

        # Combine and sort the results
        single_chrom_hmmCDR_output = pd.concat([single_chrom_CDRs, single_chrom_low_conf_CDRs])

        single_chrom_hmmCDR_output['thickStart'] = single_chrom_hmmCDR_output['start']
        single_chrom_hmmCDR_output['thickEnd'] = single_chrom_hmmCDR_output['end']

        # Assign color based on the 'name' column
        single_chrom_hmmCDR_output['itemRgb'] = np.where(
            single_chrom_hmmCDR_output['name'] == f"sub{self.output_label}", 
            self.main_color, 
            np.where(single_chrom_hmmCDR_output['name'] == f'low_conf_sub{self.output_label}', self.low_conf_color, '')
        )

        single_chrom_hmmCDR_output = single_chrom_hmmCDR_output.astype({
            'start': 'int64',
            'end': 'int64',
            'thickStart': 'int64',
            'thickEnd': 'int64'
        })

        return single_chrom_hmmCDR_output.sort_values(by='start'), single_chrom_hmmCDR_scores
    
    def hmm_single_chromosome(self, chrom, bed4Methyl_chrom, priors_chrom):
        labelled_bed4Methyl_chrom = self.assign_emissions(
            self.assign_priors(bed4Methyl_chrom, priors_chrom), 
            self.calculate_emission_thresholds(bed4Methyl_chrom)
        )

        emission_matrix = self.calculate_emission_matrix(labelled_bed4Methyl_chrom)
        transition_matrix = self.calculate_transition_matrix(labelled_bed4Methyl_chrom)

        hmmlabelled_bed4Methyl = self.runHMM(labelled_bed4Methyl_chrom, transition_matrix, emission_matrix)
        single_chrom_hmmCDR_result, single_chrom_hmmCDR_scores = self.create_subCDR_df(hmmlabelled_bed4Methyl)

        return chrom, single_chrom_hmmCDR_result, single_chrom_hmmCDR_scores

    def hmm_all_chromosomes(self, bed4Methyl_chrom_dict, priors_chrom_dict):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(self.hmm_single_chromosome, chrom, bed4Methyl_chrom_dict[chrom], priors_chrom_dict[chrom]): chrom 
                for chrom in bed4Methyl_chrom_dict
            }

            results = {chrom: future.result() for future, chrom in futures.items()}
            hmmCDRresults_chrom_dict = {chrom: result[1] for chrom, result in results.items()}
            hmmCDRscores_chrom_dict = {chrom: result[2] for chrom, result in results.items()}

        return hmmCDRresults_chrom_dict, hmmCDRscores_chrom_dict


def main():
    argparser= argparse.ArgumentParser(description='Process input files with optional parameters.')

    argparser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    argparser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    argparser.add_argument('output_path', type=str, help='Output Path for the output files')

    # hmmCDR Parser Flags
    argparser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    argparser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")')
    argparser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--min_valid_cov', type=int, default=1, help='Minimum Valid Coverage to consider a methylation site. (default: 1)')

    # hmmCDR Priors Flags
    argparser.add_argument('--window_size', type=int, default=1000, help='Window size to calculate prior regions. (default: 1000)')
    argparser.add_argument('--prior_percentile', type=float, default=5, help='Percentile for finding  prior subCDR regions. (default: 5)')
    argparser.add_argument('--prior_min_size', type=int, default=5000, help='Minimum size of prior subCDR regions. (default: 5000)')
    argparser.add_argument('--prior_merge_distance', type=int, default=1001, help='Distance to merge adjacently labelled subCDR regions. (default: 1001)')

    # HMM Flags
    argparser.add_argument('--raw_thresholds', action='store_true', default=True, help='Use values for flags w,x,y,z as raw threshold cutoffs for each emission category. (default: True)')
    argparser.add_argument('--n_iter', type=int, default=1, help='Maximum number of iteration allowed for the HMM. (default: 1)')
    argparser.add_argument('--hmm_min_size', type=int, default=1000, help='Minimum size of region identified. (default: 1000)')
    argparser.add_argument('--hmm_merge_distance', type=int, default=1001, help='Distance to merge adjacently labelled subCDR regions. (default: 1001)')
    argparser.add_argument('--min_cdr_score', type=float, default=0.95, help='The minimum HMM score [0-1] required to call a CDR. (default: 0.95)')
    argparser.add_argument('--min_low_conf_score', type=float, default=0.75, help='The minimum HMM score [0-1] required to call a low confidence CDR. (default: 0.75)')
    argparser.add_argument('--main_color', type=str, default='50,50,255', help='Color to dictate main regions. (default: 50,50,255)')
    argparser.add_argument('--low_conf_color', type=str, default='100,150,200', help='Color to dictate low confidence regions. (default: 100,150,200)')
    argparser.add_argument('-w', type=int, default=0, help='Threshold of non-zero methylation percentile to be classified as very low (default: 0)')
    argparser.add_argument('-x', type=int, default=25, help='Threshold of non-zero methylation percentile to be classified as low (default: 25)')
    argparser.add_argument('-y', type=int, default=50, help='Threshold of non-zero methylation percentile to be classified as medium (default: 50)')
    argparser.add_argument('-z', type=int, default=75, help='Threshold of non-zero methylation percentile to be classified as high (default: 75)')

    # Shared Flags
    argparser.add_argument('--enrichment', action='store_true', default=False, help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    argparser.add_argument('--merge_distance', type=int, default=300000, help='Distance to merge subCDRs into a CDR. (default: 300000)')
    argparser.add_argument('--min_subCDRs', type=int, default=3, help='Minimum number of subCDRs to report a CDR. (default: 3)')
    argparser.add_argument('--output_all', action='store_true', default=False, help="Set to true if you would like to save all intermediate filesf. (default: False)")
    argparser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "subCDR")')

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

    if args.output_all:
        concat_filtered_cenSat = pd.concat(cenSat_chrom_dict.values(), axis=0)
        concat_filtered_cenSat.to_csv(f'{output_prefix}_{args.sat_type}_regions.bed', sep='\t', index=False, header=False)
        concatenated_bed4Methyl = pd.concat(bed4Methyl_chrom_dict.values(), axis=0)
        concatenated_bed4Methyl.to_csv(f'{output_prefix}_{args.sat_type}_methylation.bedgraph', sep='\t', index=False, header=False)
    
    CDRpriors = hmmCDR_prior_finder(
        window_size=args.window_size, 
        prior_percentile=args.prior_percentile, 
        min_size=args.prior_min_size, 
        merge_distance=args.prior_merge_distance, 
        enrichment=args.enrichment, 
        output_label=args.output_label
    )

    priors_chrom_dict, prior_windows_chrom_dict = CDRpriors.priors_all_chromosomes(bed4Methyl_chrom_dict=bed4Methyl_chrom_dict)
    
    if args.output_all:
        concatenated_priors = pd.concat(priors_chrom_dict.values(), axis=0)
        concatenated_priors.to_csv(f'{output_prefix}_priorCDRs.bed', sep='\t', index=False, header=False)

        concatenated_windows = pd.concat(prior_windows_chrom_dict.values(), axis=0)
        concatenated_windows.to_csv(f'{output_prefix}_windowmeans.bedgraph', sep='\t', index=False, header=False)

    CDRhmm = hmmCDR(
        raw_thresholds=args.raw_thresholds,
        n_iter=args.n_iter,
        min_size=args.hmm_min_size,
        merge_distance=args.hmm_merge_distance, 
        min_cdr_score=args.min_cdr_score,
        min_low_conf_score=args.min_low_conf_score,
        main_color=args.main_color,
        low_conf_color=args.low_conf_color,
        w=args.w, x=args.x, y=args.y, z=args.z,
        output_label=args.output_label
    )

    hmm_results_chrom_dict, hmm_scores_chrom_dict = CDRhmm.hmm_all_chromosomes(
        bed4Methyl_chrom_dict=bed4Methyl_chrom_dict,
        priors_chrom_dict=priors_chrom_dict
    )

    # output subCDRs and CDR scoring bedgraph
    concatenated_hmm_subCDRs = pd.concat(hmm_results_chrom_dict.values(), axis=0)
    concatenated_hmm_subCDRs.to_csv(f'{output_prefix}_sub{args.output_label}.bed', sep='\t', index=False, header=False)
    concatenated_hmm_scores = pd.concat(hmm_scores_chrom_dict.values(), axis=0)
    concatenated_hmm_scores.to_csv(f'{output_prefix}_scores.bedgraph', sep='\t', index=False, header=False)

    # create final CDR output that merges adjacent subCDRs and reports scoring of over 3
    concatenated_hmm_subCDRs_bedtools = pybedtools.BedTool.from_dataframe(concatenated_hmm_subCDRs[concatenated_hmm_subCDRs.iloc[:,3] == f"sub{args.output_label}"])
    concatenated_hmm_CDRs = concatenated_hmm_subCDRs_bedtools.merge(d=args.merge_distance, c=2, o='count').to_dataframe(names=['chrom', 'start', 'end', 'count'])
    concatenated_hmm_CDRs['name'] = concatenated_hmm_CDRs.loc[:, 'count'].apply(
            lambda x: f"low_conf_{args.output_label}" if x < args.min_subCDRs else f"{args.output_label}"
        )
    
    # Format final output dataframe
    hmmCDR_names = concatenated_hmm_CDRs.loc[:,'name']
    hmmCDR_counts = concatenated_hmm_CDRs.loc[:,'count']
    concatenated_hmm_CDRs = concatenated_hmm_CDRs.iloc[:,:3]
    concatenated_hmm_CDRs['name'] = hmmCDR_names
    concatenated_hmm_CDRs['score'] = hmmCDR_counts
    concatenated_hmm_CDRs['strand'] = '.'
    concatenated_hmm_CDRs['thickStart'] = concatenated_hmm_CDRs.iloc[:,1]
    concatenated_hmm_CDRs['thickEnd'] = concatenated_hmm_CDRs.iloc[:,2]
    concatenated_hmm_CDRs['itemRgb'] = np.where(
        concatenated_hmm_CDRs['name'] == f"{args.output_label}", args.main_color, 
        np.where(concatenated_hmm_CDRs['name'] == f'low_conf_{args.output_label}', args.low_conf_color, ''))

    concatenated_hmm_CDRs.to_csv(args.output_path, sep='\t', index=False, header=False)


if __name__ == "__main__":
    main()
    