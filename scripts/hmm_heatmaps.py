import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse

def main():
    # read in the command line inputs
    parser = argparse.ArgumentParser(
        prog='hmm_heatmaps.py',
        description="""Separate candidate CDR containing reads from a bamfile""")

    parser.add_argument("-e", "--emissionMatrix",
                        required=True,
                        metavar="bed file containing estimate CDR Regions",
                        help="bed file containing estimate CDR Regions with chromosome, starting position and ending positon")
    parser.add_argument("-t", "--transitionMatrix",
                        required=True,
                        metavar="bed file containing estimate CDR Transition Regions",
                        help="bed file containing estimate CDR Transitions with chromosome, starting position and ending positon")
    parser.add_argument("-o", "--outputPrefix",
                        required=False,
                        default="hmm_heatmap",
                        metavar="output prefix",
                        help="The output prefix")
    
    args = parser.parse_args()
    emissionMatrix_path = args.emissionMatrix
    transitionMatrix_path = args.transitionMatrix
    outputPrefix = args.outputPrefix

    # read in the filepaths passed in and open them as numpy matrices
    emissionMatrix = pd.read_csv(emissionMatrix_path, 
                     sep=',', 
                     header=0, 
                     index_col=0)
    transitionMatrix = pd.read_csv(transitionMatrix_path, 
                     sep=',', 
                     header=0, 
                     index_col=0)
    

    ########################################
    ### make the emission matrix heatmap ###
    ########################################
    print( emissionMatrix )

    # Plot the heatmap with custom colors and annotations
    plt.imshow(emissionMatrix.values, cmap="Reds", vmin=0,
           vmax=1, extent=[0, emissionMatrix.shape[1], 0, emissionMatrix.shape[0]])
    
    reversed_eMatrix = np.flip(emissionMatrix.values, axis=0)
    for y in range( emissionMatrix.shape[0] ): 
        for x in range( emissionMatrix.shape[1] ): 
            plt.annotate(str(round(reversed_eMatrix[y, x], 8)), xy=(x+0.5, y+0.5), 
                        ha='center', va='center', color='black', fontsize=8)

    # Add colorbar 
    cbar = plt.colorbar(ticks=[0, 0.25, 0.5, 0.75, 1], shrink=0.5) 
    cbar.ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1'])

    # Set plot title and axis labels 
    emissionMatrixTitle = "Emission Matrix"
    plt.title(emissionMatrixTitle) 
    plt.xlabel("Hidden State") 
    plt.ylabel("Emission")
    plt.tick_params(top=False, labeltop=True, 
                    bottom=False, labelbottom=False,
                    left=False, labelleft=True)
    plt.xticks([0.5, 1.5, 2.5], ['A', 'B', 'C'])
    plt.yticks([0.5, 1.5, 2.5, 3.5], ['w', 'x', 'y', 'z'])

    # Display the plot
    emissionMatrixHeatmap = outputPrefix + '.emissionMatrixHeatmap.png'
    plt.savefig(emissionMatrixHeatmap)
    plt.clf()

    ##########################################
    ### make the transition matrix heatmap ###
    ##########################################
    print( transitionMatrix )

    # Plot the heatmap with custom colors and annotations
    plt.imshow(transitionMatrix.values, cmap="Reds", vmin=0,
           vmax=1, extent=[0, transitionMatrix.shape[1], 0, transitionMatrix.shape[0]])
    
    reversed_tMatrix = np.flip(transitionMatrix.values, axis=0)
    for i in range( transitionMatrix.shape[0] ): 
        for j in range( transitionMatrix.shape[1] ): 
            plt.annotate(str(round(reversed_tMatrix[i, j], 8)), xy=(j+0.5, i+0.5), 
                        ha='center', va='center', color='black', fontsize=8)

    # Add colorbar 
    cbar = plt.colorbar(ticks=[0, 0.25, 0.5, 0.75, 1], shrink=0.5) 
    cbar.ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1'])

    # Set plot title and axis labels 
    transitionMatrixTitle = "Transition Matrix"
    plt.title(transitionMatrixTitle) 
    plt.xlabel("From Hidden State") 
    plt.ylabel("To Hidden State")
    plt.tick_params(top=False, labeltop=True, 
                    bottom=False, labelbottom=False,
                    left=False, labelleft=True)
    plt.xticks([0.5, 1.5, 2.5], ['A', 'B', 'C'])
    plt.yticks([0.5, 1.5, 2.5], ['C', 'B', 'A'])

    # Display the plot 
    transitionMatrixHeatmap = outputPrefix + '.transitionMatrixHeatmap.png'
    plt.savefig(transitionMatrixHeatmap)
    plt.clf()

if __name__ == '__main__':
    main()
