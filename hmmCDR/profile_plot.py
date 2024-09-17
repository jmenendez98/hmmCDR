import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import pybedtools


def get_file_name(path):
    """Helper function to extract file name with extension."""
    return os.path.basename(path)

def parse_files(paths):
    paths_list = paths.split(',')
    column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCounts', 'blockSizes', 'blockStarts']
    parsed_files = {}
    for path in paths_list:
        data = pd.read_csv(path, sep='\t', header=None)
        n_cols = data.shape[1]
        data.columns = column_names[:n_cols]
        parsed_files[get_file_name(path)] = data
    return parsed_files

def find_center(cenSat, sat_names='H1L,mon'):
    if sat_names is None:
        sat_names = []
    elif isinstance(sat_names, str):
        sat_names = [name.strip() for name in sat_names.split(',')]
    mask = cenSat['name'].str.contains('|'.join(sat_names), na=False)
    cenSat_subset = cenSat[mask]
    subset_starts = cenSat_subset.groupby('chrom')['start'].min()
    subset_ends = cenSat_subset.groupby('chrom')['end'].max()
    subset_half_lengths = ((subset_ends - subset_starts) // 2)
    subset_middles = subset_starts + subset_half_lengths
    return subset_middles, subset_half_lengths

def normalize_to_center(bed, centers):
    if 'start' not in bed.columns or 'end' not in bed.columns:
        raise KeyError("'start' or 'end' column not found in the bed DataFrame.")
    bed['normalized_start'] = bed.apply(
        lambda row: row['start'] - centers.get(row['chrom'], 0), 
        axis=1
    )
    bed['normalized_end'] = bed.apply(
        lambda row: row['end'] - centers.get(row['chrom'], 0),  
        axis=1
    )
    return bed

def sep_haplotype(data, hap1_str, hap2_str):
    return data[data['chrom'].str.contains(hap1_str)], data[data['chrom'].str.contains(hap2_str)]

def create_range_bed(centers, half_size):
    """Create a BED file from the centers Series with chromosome, start, and end coordinates."""
    bed_df = pd.DataFrame({
        'chrom': centers.index,                 
        'start': centers - half_size,            
        'end': centers + half_size                
    }).astype({'start': 'int', 'end': 'int'})  
    bedtool = pybedtools.BedTool.from_dataframe(bed_df)
    return bedtool

def intersect_with_range(df, range_bedtool):
    """Intersect all dataframes in the files dictionary with range_bedtool."""
    bedtool = pybedtools.BedTool.from_dataframe(df)
    intersected_bedtool = bedtool.intersect(range_bedtool)
    intersected_df = intersected_bedtool.to_dataframe()
    return intersected_df

def smooth_bedgraph(bedgraph_df, window_size=10000):
    bedgraph_df.loc[:, 'name'] = pd.to_numeric(bedgraph_df['name'], errors='coerce')
    bedgraph_bt = pybedtools.BedTool.from_dataframe(bedgraph_df[['chrom', 'start', 'end', 'name']])
    
    smoothed_dfs = []
    for chrom, group in bedgraph_df.groupby('chrom'):
        range_start = group['start'].min()
        range_end = group['end'].max()

        windows_list = []
        for start in range(range_start, range_end, window_size):
            end = min(start + window_size, range_end)
            windows_list.append([chrom, start, end])
        windows_bt = pybedtools.BedTool(windows_list)
        
        mapped = windows_bt.map(bedgraph_bt, c=4, o='mean')
        mapped_df = mapped.to_dataframe(names=['chrom', 'start', 'end', 'name'])
        mapped_df['name'] = mapped_df['name'].replace('.', np.nan).fillna(0)
        smoothed_dfs.append(mapped_df)
    
    if smoothed_dfs:
        return pd.concat(smoothed_dfs).reset_index(drop=True)
    else:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'name'])

def cenprofileplot(hap_dict, num_features, cenSat_name, 
                   x_min, x_max, hap1_name, hap2_name,
                   output_path, no_track_labels, bedgraphs):
    haplotypes = [hap1_name, hap2_name] if hap2_name else [hap1_name]
    num_haplotypes = len(haplotypes)
    num_chroms = max([len(hap_dict[haplotype][cenSat_name]['chrom'].unique()) for haplotype in haplotypes])

    fig, axes = plt.subplots(1, num_haplotypes, figsize=(40, num_chroms*0.5*num_features))

    if num_haplotypes == 1:
        axes = [axes]

    def get_chrom_num(chrom):
        match = re.match(r'chr(\d+|X|Y)_(MATERNAL|PATERNAL)', chrom)
        if match:
            chrom_num = match.group(1)
            if chrom_num.isdigit():
                return int(chrom_num)
            elif chrom_num in ['X', 'Y']:
                return 23
        return float('inf')

    def plot_bars_and_lines(ax, bed_dictionary, num_features, y_pos_dict):
        index = -1
        for key, bed_data in bed_dictionary.items():
            unique_chroms = sorted(bed_data['chrom'].unique(), key=get_chrom_num)
            index += 1 

            for chrom in unique_chroms:
                data = bed_data[bed_data['chrom'] == chrom]
                y_pos = y_pos_dict[chrom]
                
                starts = data['normalized_start'].values
                ends = data['normalized_end'].values
                widths = ends - starts 

                if key not in bedgraphs:
                    if 'itemRgb' in data.columns:
                        colors = data['itemRgb'].apply(lambda x: tuple([int(c) / 255 for c in x.split(',')]) )
                    else:
                        colors = 'grey'
                    if key == cenSat_name:
                        bar_heights = 0.9
                        offset = 0
                    else: 
                        bar_heights = 0.35
                        offset = (index*0.4)+0.4 
                    ax.barh(
                        y_pos * num_features + offset,  
                        widths,
                        left=starts,
                        height=bar_heights,
                        color=colors
                    )
                else:
                    height_factors = (data['name'].astype(float).values / 100)
                    if key == cenSat_name:
                        bar_heights = height_factors * 0.9
                        offset = 0
                    else: 
                        bar_heights = height_factors * 0.35
                        offset = (index*0.4)+0.4 
                    ax.barh(
                        y_pos * num_features + offset, 
                        widths,
                        left=starts,
                        height=bar_heights,
                        color='grey',
                        align='edge'
                    )

                if not no_track_labels:
                    ax.text(
                        x_min + 0.001 * (x_max - x_min), 
                        y_pos * num_features + offset,
                        key, 
                        va='center',
                        ha='left',   
                        fontsize=10, 
                        fontweight=1000,
                        color='black'
                    )

    for ax, haplotype in zip(axes, haplotypes):
        hap_cenSat = hap_dict[haplotype][cenSat_name]
        hap_chrom_labels = sorted(hap_cenSat['chrom'].unique(), key=get_chrom_num)[::-1]
        y_pos_dict = {chrom: i for i, chrom in enumerate(hap_chrom_labels)}

        plot_bars_and_lines(ax, hap_dict[haplotype], num_features, y_pos_dict)

        ax.set_yticks([i * num_features + num_features / 4 for i in range(len(hap_chrom_labels))])
        ax.set_yticklabels(hap_chrom_labels, fontsize=14)

        ax.set_ylim(-0.5, len(hap_chrom_labels) * num_features - 0.5)
        ax.set_xlim(x_min, x_max)
        ax.xaxis.set_ticks([])

    plt.tight_layout()
    plt.savefig(output_path, dpi=600)

def main():

    def parse_arguments():
        parser = argparse.ArgumentParser(description="Process genomic data.")

        parser.add_argument('cenSat_path', type=str,
                            help='Path to the cenSat file.')
        parser.add_argument('output_path', type=str, 
                            help='Path to the bedgraph file.')
        
        parser.add_argument('--bed_paths', type=str, default='',
                            help='Comma-separated list of paths to the bed files.')
        parser.add_argument('--bedgraph_paths', type=str, default='',
                            help='Path to the bedgraph file.')
        parser.add_argument('--hap1_name', type=str, default='MATERNAL',
                            help='Name of the first haplotype.')
        parser.add_argument('--hap2_name', type=str, default='PATERNAL',
                            help='Name of the second haplotype.')
        parser.add_argument('--diploid', action='store_true', default=True,
                            help='Set to True if the data is diploid.')
        parser.add_argument('--no_track_labels', action='store_true', default=False,
                    help='Set to True if the data is diploid.')

        return parser.parse_args()

    args = parse_arguments()

    cenSat_path = args.cenSat_path
    bed_paths = args.bed_paths
    bedgraph_paths = args.bedgraph_paths
    output_path = args.output_path
    hap1_name = args.hap1_name
    hap2_name = args.hap2_name
    diploid = args.diploid
    no_track_labels = args.no_track_labels

    cenSat_name = get_file_name(cenSat_path)

    paths = ','.join([cenSat_path, bed_paths, bedgraph_paths])

    files = parse_files(paths)

    centers, half = find_center(files[cenSat_name], sat_names='H1L')
    x_min, x_max = -(half.max()*1.1), (half.max()*1.1)

    range_bedtool = create_range_bed(centers, half_size=(half.max()*1.1) )

    files = {file: intersect_with_range(files[file], range_bedtool) for file in files}

    if diploid:
        hap_dict = { hap1_name: {}, hap2_name: {} }

        for file_name in files.keys():
            hap1_file, hap2_file = sep_haplotype(files[file_name], hap1_name, hap2_name)
            if not file_name.endswith('.bedgraph'):
                hap_dict[hap1_name][file_name], hap_dict[hap2_name][file_name] = hap1_file, hap2_file
            else:
                hap_dict[hap1_name][file_name] = smooth_bedgraph(hap1_file, window_size=2500) 
                hap_dict[hap2_name][file_name] = smooth_bedgraph(hap2_file, window_size=2500) 
                
    else:
        hap1_name = 'haploid'
        hap_dict = { hap1_name: {} }

        for file_name in files.keys():
            if file_name.endswith('.bedgraph'):
                hap_dict[hap1_name][file_name] = files[file_name]
            else:
                hap_dict[hap1_name][file_name] = smooth_bedgraph(files[file_name], window_size=2500) 

    for haplotype in hap_dict.keys():
        for bed in hap_dict[haplotype].keys():
            hap_dict[haplotype][bed] = normalize_to_center(hap_dict[haplotype][bed], centers=centers)

    cenprofileplot(hap_dict=hap_dict, num_features=len(files.keys()), 
               cenSat_name=cenSat_name, hap1_name=hap1_name, hap2_name=hap2_name,
               x_min=x_min, x_max=x_max,
               output_path=output_path, no_track_labels=no_track_labels, 
               bedgraphs=bedgraph_paths)
    

if __name__ == "__main__":
    main()