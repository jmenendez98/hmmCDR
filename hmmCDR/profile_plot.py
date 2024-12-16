import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import re
import pybedtools

def get_file_name(path):
    """Helper function to extract file name with extension."""
    return os.path.basename(path)

def parse_files(paths):
    """Parse multiple file paths into DataFrames."""
    paths_list = paths.split(',')
    return {get_file_name(path): pybedtools.BedTool(path) for path in paths_list if path}

def find_center(censat_bedtool, sat_names='active_hor'):
    """Find chromosome centers using BedTool operations."""
    # Convert to dataframe only once for the groupby operation
    censat = censat_bedtool.to_dataframe()
    
    # Handle sat_names as either string or list
    if isinstance(sat_names, str):
        sat_names = [name.strip() for name in sat_names.split(',')]
    
    # Use proper pandas string method and column name
    # BedTool columns are named 'name' for the 4th column, not 3
    mask = censat['name'].str.contains('|'.join(sat_names), na=False)
    censat = censat[mask]

    subset_starts = censat.groupby('chrom')['start'].min()
    subset_ends = censat.groupby('chrom')['end'].max()
    subset_half_lengths = ((subset_ends - subset_starts) // 2)
    return subset_starts + subset_half_lengths, subset_half_lengths

def normalize_to_center(bed, centers):
    """Normalize bed coordinates relative to chromosome centers."""
    bed = bed.copy()
    bed['normalized_start'] = bed['start'] - bed['chrom'].map(centers).fillna(0)
    bed['normalized_end'] = bed['end'] - bed['chrom'].map(centers).fillna(0)
    return bed

def create_range_bed(centers, half_size):
    """Create a BED file from chromosome centers."""
    bed_df = pd.DataFrame({
        'chrom': centers.index,                 
        'start': centers - half_size,            
        'end': centers + half_size                
    }).astype({'start': 'int', 'end': 'int'})  
    return pybedtools.BedTool.from_dataframe(bed_df)

def cen_profile_plot(files, output_file, censat_name, centers_dict, x_min, x_max, no_track_labels=False):
    """
    Create summary centromeric profile plot across genome.
    
    Parameters:
    -----------
    files : list
        List of input files
    censat_name : str
        Name of the centromere satellite file
    x_min : float
        Minimum x-axis value
    x_max : float
        Maximum x-axis value
    no_track_labels : bool, optional
        Whether to show track labels (default: False)
    """
    
    plt.figure(figsize=(12,24))

    def get_chrom_num(chrom):
        """Extract chromosome number for sorting."""
        match = re.match(r'chr(\d+|X|Y)(_MATERNAL|_PATERNAL)?', chrom)
        if match:
            chrom_num = match.group(1)
            if chrom_num.isdigit():
                return int(chrom_num)
            elif chrom_num == 'X':
                return 23
            elif chrom_num == 'Y':
                return 24
        return float('inf')

    all_chromosomes = sorted(
        set(sorted(chrom for file in files.values() for chrom in file.to_dataframe()['chrom'].unique())),
        key=get_chrom_num
    )

    plt.xlim(x_min, x_max)
    plt.ylim(0, len(all_chromosomes)*len(files)+1)
    ax = plt.gca()    

    def parse_rgb_color(rgb_str):
        """Convert RGB string to color tuple."""
        try:
            # Handle both comma-separated and space-separated RGB values
            rgb = [int(x)/255.0 for x in rgb_str.replace(',', ' ').split()]
            return tuple(rgb[:3] + [0.7])  # Add alpha for consistency
        except (ValueError, IndexError):
            return None

    def plot_bed(file, file_name, chrom_order, centers_dict, offset, no_track_label=no_track_labels):
        """
        Plot genomic data tracks from multiple BED files.
        
        Parameters:
        -----------
        file : pybedtools BedTool
            BED file data
        chrom_order : list
            Sorted list of chromosomes
        no_track_label : bool, optional
            Whether to show track labels
        """
        data = file.to_dataframe()

        for chrom in chrom_order:
            # Find the index of the chromosome to determine y-coordinate
            # Subtract from len(chrom_order) to have chr1 at the top
            y = ( len(chrom_order) - chrom_order.index(chrom) ) * len(files) - len(files)+1 + offset

            if not no_track_label:
                plt.text(x_min, y, file_name, va='bottom', ha='left', size=4)
            
            # Filter data for this specific chromosome
            chrom_data = data[data['chrom'] == chrom]
            
            # Get the centromere position for this chromosome
            center = centers_dict[chrom]

            for i, row in chrom_data.iterrows():
                start = int(row.iloc[1]) - int(center)
                length = int(row.iloc[2]) - int(row.iloc[1])
                try:
                    color = parse_rgb_color(row.iloc[8])
                except IndexError:
                    color = 'black'
                
                rect = patches.Rectangle((start, y), length, 0.5, color=color)
                ax.add_patch(rect)

    offset=0
    for file in files:
        plot_bed(files[file], file, all_chromosomes, centers_dict, offset)
        offset += 1

    # Create y-tick labels for chromosomes
    y_ticks = []
    y_tick_labels = []
    for i, chrom in enumerate(all_chromosomes):
        # Calculate y position for each chromosome
        y_pos = (len(all_chromosomes) - i) * len(files) - len(files)//2
        y_ticks.append(y_pos)
        y_tick_labels.append(chrom)

    plt.yticks(y_ticks, y_tick_labels, fontsize=8)

    plt.savefig(output_file, dpi=1200)


def main():
    # argparse files in...
    argparser= argparse.ArgumentParser(description='Process input files with optional parameters.')

    argparser.add_argument('censat', type=str, help='Path to the bedMethyl file')
    argparser.add_argument('output', type=str, help='Path (with desired extension) of output file')
    argparser.add_argument('-b', '--beds', type=str, default='', help='Paths to the BED files (comma separated).')
    argparser.add_argument('-g', '--bedgraphs', type=str, default='', help='Paths to the BEDGRAPH files (comma separated). [WIP]')

    args = argparser.parse_args()

    censat_name = get_file_name(args.censat)
    bedgraph_names = get_file_name(args.bedgraphs)

    paths = ','.join([args.censat, args.beds, args.bedgraphs])
    files = parse_files(paths)

    centers, half = find_center(files[censat_name], sat_names='H1L')
    x_min, x_max = -(half.max()*1.1), (half.max()*1.1)

    range_bedtool = create_range_bed(centers, half_size=(half.max()*1.1) )
    files = {file: files[file].intersect(range_bedtool, wa=True) for file in files}

    cen_profile_plot(files=files, output_file=args.output, censat_name=censat_name, centers_dict=centers, x_min=x_min, x_max=x_max)


if __name__ == "__main__":
    main()