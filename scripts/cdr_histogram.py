import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create histogram of modification probability distribution with colored bars representing region proportions.")
    parser.add_argument("-i", "--modification_probabilities", 
                        help="Path to the modification probabilities bed file",
                        required=True)
    parser.add_argument("-r", "--regions", 
                        help="Path to the regions bed file",
                        required=True)
    parser.add_argument("-o", "--output", 
                        help="Path to save the histogram plot", 
                        default="histogram.png")
    return parser.parse_args()

def calculate_bins(mod_prob_df, regions_df, hist_edges):
    # Initialize list to store counts for each bin
    bin_counts = []
    bin_proportions = []

    region_types = ['Not within region', 'CDR', 'CDR_Transition', 'small_CDR']

    # Iterate over histogram bins
    for i in range(len(hist_edges) - 1):
        start, end = hist_edges[i], hist_edges[i + 1]
        
        # Filter entries within the current bin's range
        entries_in_bin = mod_prob_df[(mod_prob_df[3] >= start) & (mod_prob_df[3] < end)]

        # Count the number of entries in the current bin and append to bin_counts
        bin_counts.append(len(entries_in_bin))

        # Initialize a dictionary to store counts for each region type
        region_counts = {region_type: 0 for region_type in region_types}

        for _, row in entries_in_bin.iterrows():
            placed = False
            #print("ROW:", row)
            for _, region in regions_df.iterrows():
                #print("REGION:", region)

                # Check if the entry falls within a region and update counts accordingly
                if region[1] <= row[1] < region[2]:
                    region_counts[region[3]] += 1  # Increment count for the region type
                    placed = True
                    break
            
            if not placed:
                region_counts['Not within region'] += 1

        bin_proportions.append(region_counts)

    return bin_counts, bin_proportions

def main():
    args = parse_arguments()

    # Read bed files into pandas DataFrames
    mod_prob_df = pd.read_csv(args.modification_probabilities, sep="\t", header=None)
    regions_df = pd.read_csv(args.regions, sep="\t", header=None)

    # Generate histogram edges from 0 to 100 with step size 5
    hist_edges = np.arange(0, 101, 5)

    # Calculate histogram using the specified edges
    hist, _ = np.histogram(mod_prob_df[3], bins=hist_edges, density=True)

    #print("Histogram Edges:", hist_edges)
    #print("ModProb DF:", mod_prob_df)
    #print("Regions DF:", regions_df)

    # Calculate region proportions for each histogram box
    bin_counts, bin_proportions = calculate_bins(mod_prob_df, regions_df, hist_edges)

    #print("Bin Counts:", bin_counts)
    #print("Bin Proportions:", bin_proportions)

    # Define colors for each region
    region_colors = {'Not within region': 'lightsteelblue', 
                     'CDR': 'maroon', 
                     'CDR_Transition': 'wheat', 
                     'small_CDR': 'orangered'}

    #print("Bin Counts:", bin_counts)
    #print("Bin Proportions:", bin_proportions)

    # Define region order
    region_order = ['CDR', 
                    'CDR_Transition',
                    'small_CDR', 
                    'Not within region']

    # Convert bin counts and proportions to numpy arrays
    bin_counts = np.array(bin_counts)
    bin_proportions = np.array([list(proportions.values()) for proportions in bin_proportions])

    # Create a vertical bar plot
    fig, ax = plt.subplots()

    # Initialize the bottom of each bar to be plotted
    bottom = np.zeros_like(bin_counts)

    index = {'Not within region': 0, 
             'CDR': 1, 
             'CDR_Transition': 2, 
             'small_CDR': 3}

    # Iterate over each region type
    for region_type in region_order:
        #print(region_type, index[region_type])
        # Plot the stacked bars for the current region type
        ax.bar(x=hist_edges[:-1], 
               height=bin_proportions[:, index[region_type]], 
               width=(5),
               bottom=bottom, 
               color=region_colors[region_type], 
               label=region_type,
               align='edge')
        # Update the bottom of the next region to be plotted
        bottom += bin_proportions[:, index[region_type]]

    # Set y-axis label
    ax.set_ylabel('Counts')

    # Set y-axis limits based on the total number of regions in each bin
    ax.set_ylim(0, bin_proportions.sum(axis=1).max() * 1.1)

    # Set x-axis label
    ax.set_xlabel('Modification Probability')

    # Set x-axis limits
    ax.set_xlim(0, 100)

    # Set x-ticks
    ax.set_xticks(hist_edges)

    # Set only every fourth tick label as visible
    for label in ax.xaxis.get_ticklabels():
        tick_index = int(float(label.get_text()))
        if tick_index % 10 == 0 or tick_index == 0:
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Set title
    ax.set_title('Histogram of Region Counts per Bin')

    # Add legend
    ax.legend()

    # Save plot as an image
    plt.savefig(args.output, dpi=3200)


if __name__ == "__main__":
    main()

