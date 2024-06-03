import argparse
import numpy as np
import pandas as pd

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        self.parser = argparse.ArgumentParser(
            prog='HMMCDRReferenceDetection.py',
            description="""Identify High Confidence CDR Regions and Transitions Using an HMM""")
        
        self.parser.add_argument("-l", "--learningRate",
                            action = 'store', 
                            nargs='?', 
                            const=True,
                            required=False,
                            default="0.0000001",
                            metavar="Learning Rate for the Viterbi Learning",
                            help="Sets the learning rate for the Viterbi Learning")
        self.parser.add_argument("--maxSteps",
                            action = 'store', 
                            nargs='?', 
                            const=True,
                            required=False,
                            default="100",
                            metavar="Maximum number of Viterbi Learning Iterations",
                            help="Sets the maximum number of Viterbi Learning Iterations")
        self.parser.add_argument("-p", "--modPosFile",
                            required=True,
                            metavar="bed4 file containing modified CpG site fraction modified",
                            help="bed4 file containing modified CpG site fraction modified")
        self.parser.add_argument("-s", "--strictCDRFile",
                            required=True,
                            metavar="bed File containing Initial CDR Site Estimates(from strictCDR)",
                            help="bed File containing Initial CDR Site Estimates(from strictCDR)")
        self.parser.add_argument("-o", "--outputPrefix",
                            required=True,
                            metavar="output bed file",
                            help="The output bed file name")
        
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class PathEstimate:
    '''
    Creates the initial path estimate from CG positions from reference genome for Viterbi Leanring. Also
    saves other important variable such as CpG positions and their methylation probabilities and 
    current chromosome

    Attributes:
    - cpgSitesHORFile: Subsetted bed file containing CpG sites and probabilties in HOR regions 
    - cdrStrictRegionsFile: Bed file rough estimates of CDR regions

    Parses data in these files to create initial estimates of transition and emission probabilities
    '''
    def __init__ (self, cpgSitesHORFile):
        '''
        Store HOR CpG site probability bed files
        '''
        self.self.cpgSitesHORFile = cpgSitesHORFile
    def getPath (self):
        '''
        Creates the initial path estimate for initial Viterbi Learning

        Addtionally creates global variables such as:
        - self.chrName: Chromosome label/name from bed files
        - self.cpgSitesAndProbs: Dictionary with key as CpG starting postion and value as methylation probability
        '''
        # Creates variables to store CG postions, chromosome names, and estimated path
        self.cpgSitesAndProbs = {}
        self.chrName = ""
        path = ""
        # Opens the CpG HOR methylation probability bed file as input file
        with open(self.cpgSitesHORFile) as inputFile:
            # Loops through every line in the input file
            for line in inputFile:
                # Splits line into columns in a list
                fileColumnSplit = line.strip().split("\t")
                # Saves chromosome name in column 1 of the file
                self.chrName = fileColumnSplit[0]
                # Saves CpG site starting position and methylation probability
                self.cpgSitesAndProbs[int(fileColumnSplit[1])] = float(fileColumnSplit[3])
                # Checks if CpG methylation probability is greater than 50% to store as methylation count
                if float(fileColumnSplit[3]) >= 50:
                    path += "x"
                else:
                    path += "y"
        return path

    def getChrName(self):
        '''
        Returns chromosome name
        '''
        return self.chrName
    
    def getCpGSitesAndProbs(self):
        '''
        Returns dictionary of CpG starting positions as keys and methylation probabilities as values
        '''
        return self.cpgSitesAndProbs



class InitialMatricesEstimate:
    '''
    Creates the initial transition and emission matrices for Viterbi Learning of CDR Regions. Also
    saves other important variable such as CpG positions and their methylation probabilities and 
    current chromosome

    Attributes:
    - cpgSitesHORFile: Subsetted bed file containing CpG sites and probabilties in HOR regions 
    - cdrStrictRegionsFile: Bed file rough estimates of CDR regions

    Parses data in these files to create initial estimates of transition and emission probabilities
    '''
    def __init__ (self, cpgSitesHORFile, cdrStrictRegionsFile):
        '''
        Store CDR region and HOR CpG site probability bed files
        '''
        self.cpgSitesHORFile = cpgSitesHORFile
        self.cdrStrictRegionsFile = cdrStrictRegionsFile

    def getInitialTransitionMatrix (self):
        '''
        Creates the initial transition matrix for CDR Viterbi Learning

        Addtionally creates global variables such as:
        - self.cdrRegions: List of tuples where each tuple represents a CDR region
        - self.cdrTransitions: List of tuples where each tuple represents a CDR transition region
        - self.chrName: Chromosome label/name from bed files
        - self.cpgSitesAndProbs: Dictionary with key as CpG starting postion and value as methylation probability
        '''

        # Creates variables to store CDR region postions and chromosome names
        self.cdrRegions = []
        self.cdrTransitions = []
        self.chrName = ""

        # Opens the Strict CDR File as input file
        with open(self.cdrStrictRegionsFile) as inputFile:
            # Loops through each line in the file
            for line in inputFile:
                # Splits line into columns in a list
                fileColumnSplit = line.strip().split("\t")
                self.chrName = fileColumnSplit[0]

                if fileColumnSplit[3] == 'strict_CDR':
                    # Saves CDR Regions from starting and ending positions
                    self.cdrRegions.append((int(fileColumnSplit[1]), int(fileColumnSplit[2])))
                elif fileColumnSplit[3] == 'strict_Transition':
                    # Saves CDR Regions from starting and ending positions
                    self.cdrTransitions.append((int(fileColumnSplit[1]), int(fileColumnSplit[2])))

        # Creates variables to store the initial transition matrix, counts of CpG sites in 
        # CDR/non-CDR positons and CpG sites starting positions along with methylation probabilities
        transitionMatrix = {'AA': 1e-10, 'AB': 1e-10, 'AC': 1e-10, 
                            'BA': 1e-10, 'BB': 1e-10, 'BC': 1e-10,
                            'CA': 1e-10, 'CB': 1e-10, 'CC': 1e-10,}
        prevStateCounts = {'A': 2e-10, 'B': 2e-10, 'C': 2e-10}
        self.cpgSitesAndProbs = {}

        # Opens the CpG HOR methylation probability bed file as input file
        with open(self.cpgSitesHORFile) as inputFile:
            # Loops through every line in the input file
            for line in inputFile:
                # Splits line into columns in a list
                fileColumnSplit = line.strip().split("\t")
                # Saves CpG site starting position and methylation probability
                try:
                    self.cpgSitesAndProbs[int(fileColumnSplit[1])] = float(fileColumnSplit[3])
                except:
                    continue

        # Set variable thresholds for the different types of emissions (2/4)
        # first value is 0, next is Q33 of non-zeros, last is Q66 of non-zeros (zeros are the majority most of the time)
        nonzeros = sorted( [value for value in self.cpgSitesAndProbs.values() if value != 0] )
        self.cpgThresholds = [0.0,
                              np.percentile( nonzeros, 100/3 ),
                              np.percentile( nonzeros, 200/3 )]

        # Variables to loop through CDR regions based on CpG site position and previous CpG site
        # in a CDR or not in a CDR
        cdrIndex = 0
        transitionIndex = 0
        prevCPGsiteCDR = None

        # Loops through CpG sites
        for posNum, cpgSitePos in enumerate(self.cpgSitesAndProbs.keys()):
            state = 'A'  # Assume non-CDR ('A') by default

            # Check for CDR region
            while cdrIndex < len(self.cdrRegions) and cpgSitePos > self.cdrRegions[cdrIndex][1]:
                cdrIndex += 1
            if cdrIndex < len(self.cdrRegions) and self.cdrRegions[cdrIndex][0] <= cpgSitePos <= self.cdrRegions[cdrIndex][1]:
                state = 'C'  # In CDR region
            
            # Check for CDR transition region
            while transitionIndex < len(self.cdrTransitions) and cpgSitePos > self.cdrTransitions[transitionIndex][1]:
                transitionIndex += 1
            if transitionIndex < len(self.cdrTransitions) and self.cdrTransitions[transitionIndex][0] <= cpgSitePos <= self.cdrTransitions[transitionIndex][1]:
                state = 'B'  # In CDR transition region

            # Update transition matrix
            if posNum != 0:
                prevState = prevCPGsiteState
                transitionMatrix[f'{prevState}{state}'] += 1
                prevStateCounts[prevState] += 1

            prevCPGsiteState = state

        # Correct emission matrix to have all CDRs enter/exit to a transition first(no matter how small)
        transitionMatrix['BA'] = transitionMatrix['BA'] + transitionMatrix['CA']
        transitionMatrix['CA'] = 1e-20
        transitionMatrix['BC'] = transitionMatrix['BC'] + transitionMatrix['AC']
        transitionMatrix['AC'] = 1e-20

        # Normalize transition matrix
        for key in transitionMatrix:
            prevState = key[0]
            transitionMatrix[key] /= prevStateCounts[prevState]
        
        #print('transitionMatrix:', transitionMatrix)
        return transitionMatrix
    
    def getInitialEmissionMatrix (self):
        '''
        Creates the initial emission matrix for CDR Viterbi Learning
        '''

        # Stores emission matrix, # of CpG sites in CDR/not in CDRs, 
        # path/methylation across CpG sites and index of CDR examined for overlap
        emissionMatrix = {'Aw': 1e-10, 'Ax': 1e-10, 'Ay': 1e-10, 'Az': 1e-10, 
                          'Bw': 1e-10, 'Bx': 1e-10, 'By': 1e-10, 'Bz': 1e-10, 
                          'Cw': 1e-10, 'Cx': 1e-10, 'Cy': 1e-10, 'Cz': 1e-10,}
        emissionStateCounts = {'A': 2e-10, 'B': 2e-10, 'C': 2e-10}
        path = ""
        cdrIndex = 0
        transitionIndex = 0

        # Loop through CpG sites
        for cpgSitePos, cpgProb in self.cpgSitesAndProbs.items():
            state = 'A'  # Assume non-CDR ('B') by default

            # Check for CDR region
            while cdrIndex < len(self.cdrRegions) and cpgSitePos > self.cdrRegions[cdrIndex][1]:
                cdrIndex += 1
            if cdrIndex < len(self.cdrRegions) and self.cdrRegions[cdrIndex][0] <= cpgSitePos <= self.cdrRegions[cdrIndex][1]:
                state = 'C'  # In CDR region

            # Check for CDR transition region
            while transitionIndex < len(self.cdrTransitions) and cpgSitePos > self.cdrTransitions[transitionIndex][1]:
                transitionIndex += 1
            if transitionIndex < len(self.cdrTransitions) and self.cdrTransitions[transitionIndex][0] <= cpgSitePos <= self.cdrTransitions[transitionIndex][1]:
                state = 'B'  # In CDR transition region

            # Update emission matrix based on methylation probability
                
            methylationState = 'w' if cpgProb <= self.cpgThresholds[0] else ('x' if cpgProb <= self.cpgThresholds[1] else ('y' if cpgProb <= self.cpgThresholds[2] else 'z') )
            emissionMatrix[state + methylationState] += 1
            emissionStateCounts[state] += 1
            path += methylationState

        # Normalize emission counts to probabilities
        for emission, emissionCount in emissionMatrix.items():
            emissionMatrix[emission] /= emissionStateCounts.get(emission[0])

        #print('transitionMatrix:', emissionMatrix)
        return path, emissionMatrix
    
    def getChrName(self):
        '''
        Returns chromosome name
        '''
        return self.chrName
    
    def getCpGSitesAndProbs(self):
        '''
        Returns dictionary of CpG starting positions as keys and methylation probabilities as values
        '''
        return self.cpgSitesAndProbs

class ViterbiLearning:
    '''
    Performs Viterbi Learning given an emission path/CpG methylation, transition matrix, and emission
    matrix to return the most probable CDR regions

    Attributes :
    - data: list of various attributes
        - data[0]: string of emission path for CpG sites (X - methylation, Y - not a methylation)
        - data[1]: list of emission states (X - methylation, Y - not a methylation)
        - data[2]: list of transition states (A - CDR Region, B - Not a CDR Region)
        - data[3]: dictionary of transition states and their probabilities
        - data[4]: dictionary of emission states and their probabilities
    '''
    def __init__(self, data):
        '''
        Store path, emission states, transition states, transition matrix and emission matrix as global variables
        '''
        self.path = data[0]
        self.emissionStates = data[1]
        self.transitionStates = data[2]
        self.transitionMatrix = data[3]
        self.emissionMatrix = data[4]
        self.learningRate = data[5]
        self.max_steps = data[6]

    def parameterEstimation(self, hiddenStates, learningRate):
        '''
        Performs parameter estimation given an emission path and hidden state path. Recalculates
        the transition matrix and emission matrix from the current emission path and hidden state path estimates.

        Attributes:
        - hiddenStates: String representing the current hidden state path estimate.
        '''
        emission_learning_rate = learningRate
        transition_learning_rate = learningRate * 0.1

        # Initialize new transition and emission matrices with small non-zero values to avoid division by zero
        newTransitionMatrix = {transition: 1e-10 for transition in self.transitionMatrix.keys()}
        newEmissionMatrix = {emission: 1e-10 for emission in self.emissionMatrix.keys()}
        # Initialize counters for transitions and emissions for normalization
        transitionTotalCount = {state: 2e-10 for state in self.transitionStates}
        emissionTotalCount = {state: 2e-10 for state in self.transitionStates}

        # Loop through the emission path to update transition and emission counts
        for currIndex in range(len(self.path) - 1):
            # Construct transition and emission strings for the current step
            newTransition = hiddenStates[currIndex:currIndex + 2]
            newEmission = hiddenStates[currIndex] + self.path[currIndex]

            # Increment counts for the observed transition and emission
            newTransitionMatrix[newTransition] += 1
            newEmissionMatrix[newEmission] += 1
            # Increment total counts for normalization
            transitionTotalCount[hiddenStates[currIndex]] += 1
            emissionTotalCount[hiddenStates[currIndex]] += 1

        # Handle the last emission state
        newEmission = hiddenStates[-1] + self.path[-1]
        newEmissionMatrix[newEmission] += 1
        emissionTotalCount[hiddenStates[-1]] += 1

        # Normalize transition probabilities
        for transition, count in newTransitionMatrix.items():
            startState = transition[0]  # Extract the start state from the transition string
            newTransitionMatrix[transition] /= transitionTotalCount[startState]
            newTransitionMatrix[transition] = ((1 - transition_learning_rate) * self.transitionMatrix[transition]) + (transition_learning_rate * newTransitionMatrix[transition])

        # Normalize emission probabilities
        for emission, count in newEmissionMatrix.items():
            state = emission[0]  # Extract the state from the emission string
            newEmissionMatrix[emission] /= emissionTotalCount[state]
            newEmissionMatrix[emission] = ((1 - emission_learning_rate) * self.emissionMatrix[emission]) + (emission_learning_rate * newEmissionMatrix[emission])


        # Return the updated transition and emission matrices
        return newTransitionMatrix, newEmissionMatrix

    def estimateBestPath (self):
        '''
        Calculates the best/most probable hidden state path from current transition and emission matrices
        using the Viterbi Algorithm. Essentially calculates the most probable path at each CpG site for 
        each ending state (A - CDR, B - Not a CDR, C - CDR Transtion)
        '''
        # Initialize a list to store the probabilities of the most probable paths at each step
        viterbiGraph = []

        # Loop through each position in the emission path (each CpG site)
        for emissionStep in range(len(self.path)):
            # Dictionary to store the most probable path and its probability for the current step
            currStepProbs = {}

            # Iterate over each possible hidden state (A, B, C)
            for indivState in self.transitionStates:  # Ensure 'C' is included in self.transitionStates
                # Check if we are at the first position in the emission path
                if emissionStep == 0:
                    # Initialize the probability for starting in each state, assuming uniform distribution
                    initProb = np.log(1.0 / len(self.transitionStates))  # Adjust if initial probabilities are not uniform
                    # Construct the emission string for the current state and observed emission
                    emissionStr = indivState + self.path[emissionStep]
                    # Set the initial probability for the current state, including the emission probability
                    currStepProbs[indivState] = initProb + np.log(self.emissionMatrix.get(emissionStr, -np.inf))  # Use a small log prob for missing entries
                else:
                    # Variables to track the best path and its probability up to the current step
                    bestCurrPath, bestProb = "", float("-inf")

                    # Iterate over paths and their probabilities from the previous step
                    for prevPath, prevProb in viterbiGraph[0].items():
                        # Construct the transition string from the last state of the previous path to the current state
                        transmissionStr = prevPath[-1] + indivState
                        # Construct the emission string for the current state and observed emission
                        emissionStr = indivState + self.path[emissionStep]
                        # Calculate the total probability for the current path candidate
                        currTotProb = prevProb + np.log(self.transitionMatrix.get(transmissionStr, -np.inf)) + np.log(self.emissionMatrix.get(emissionStr, -np.inf))
                        
                        # Update the best path and its probability if the current path candidate has a higher probability
                        if currTotProb > bestProb:
                            bestProb = currTotProb
                            bestCurrPath = prevPath + indivState
                    
                    # Store the most probable path and its probability for the current step
                    currStepProbs[bestCurrPath] = bestProb
            
            # Add the current step's probabilities to the Viterbi graph
            viterbiGraph = []
            viterbiGraph.append(currStepProbs)

        # After processing all positions, return the path with the highest probability at the last step
        return max(viterbiGraph[-1], key=viterbiGraph[-1].get)

    def performViterbiLearning(self):
        '''
        Perform Viterbi Learning until convergence of transition and emission matrices. Involves a loop
        to use emission path and transition & emission matrices to estimated best hidden state path and
        then emission path and the estimated hidden path to calculate new transition & emission matrices
        '''

        # Creates a set of transition and emission probabilties to compare if they are changing
        prevTransitionSum = set(self.transitionMatrix.values())
        prevEmissionSum = set(self.emissionMatrix.values())
        currTransitionSum = set()
        currEmissionSum = set()
        
        # Loops until the transition and emission matrices stay the same/converge
        learnCount = 0
        while prevTransitionSum != currTransitionSum or prevEmissionSum != currEmissionSum:
            prevTransitionSum = set(self.transitionMatrix.values())
            prevEmissionSum = set(self.emissionMatrix.values())

            hiddenStates = self.estimateBestPath() # Viterbi Algorithm
            self.transitionMatrix, self.emissionMatrix = self.parameterEstimation(hiddenStates, self.learningRate) # Estimate HMM Parameters

            currTransitionSum = set(self.transitionMatrix.values())
            currEmissionSum = set(self.emissionMatrix.values())

            if learnCount+1 >= self.max_steps:
                break
            learnCount += 1

        return self.estimateBestPath()

    def generateBedFile(self, chr, estimatedStates, cpgSitesAndProbs, outputPrefix):
        '''
        Creates bed file that contains CDR regions found through Viterbi Learning

        Attributes:
        - chr: string of chromosome name/label
        - estimatedStates: string of final hidden state path found in Viterbi Learning
        - cpgSitesAndProbs: dictionary of CpG starting positions as keys and their methylation probabilties as values
        '''

        newCDRRegions = [] # List to store starting and stopping positions of Viterbi Learning CDRs
        newCDRTransitions = []  # List to store CDR transition regions

        # Stores if previous and current CpG sites were in a CDR or not in a CDR
        prevState = None
        prevPos = None

        # Loop through hidden state path and starting positions of CpG sites
        for stateEst, cpgPos in zip(estimatedStates, cpgSitesAndProbs.keys()):
            currState = stateEst

            # Check for starting position of a CDR
            if prevState != "C" and currState == "C":
                newCDRRegions.append([chr, cpgPos, -1])
            
            # Check for ending position of a CDR
            elif prevState == "C" and currState != "C":
                if newCDRRegions:
                    newCDRRegions[-1][2] = prevPos  # Update the end position of the last CDR
            
            # Check for starting position of a CDR Transition
            if prevState != "B" and currState == "B":
                newCDRTransitions.append([chr, cpgPos, -1])
            
            # Check for ending position of a CDR Transition
            elif prevState == "B" and currState != "B":
                if newCDRTransitions:
                    newCDRTransitions[-1][2] = prevPos  # Update the end position of the last CDR Transition

            prevState = currState
            prevPos = cpgPos

        # Ensure the last region is closed properly
        if newCDRRegions and newCDRRegions[-1][2] == -1:
            newCDRRegions[-1][2] = list(cpgSitesAndProbs.keys())[-1]
        if newCDRTransitions and newCDRTransitions[-1][2] == -1:
            newCDRTransitions[-1][2] = list(cpgSitesAndProbs.keys())[-1]
            
        # Filter Transitions that are not CDR adjacent
        def filter_transitions(transitions, cdrs):
            filtered_transitions = []
            for chrom_t, start_t, end_t in transitions:
                for chrom_cdr, start_cdr, end_cdr in cdrs:
                    if chrom_t == chrom_cdr:
                        if (abs(start_t - start_cdr) <= 1000) or (abs(end_t - end_cdr) <= 1000) or (abs(end_t - start_cdr) <= 1000) or (abs(start_t - end_cdr) <= 1000):
                            filtered_transitions.append([chrom_t, start_t, end_t])
                            break  # Break the inner loop if transition is within 1000 bp of a CDR
            return filtered_transitions
        newCDRTransitions = filter_transitions(transitions=newCDRTransitions, cdrs=newCDRRegions)

        # Output CDR regions to a BED file
        output_lines = []
        for cdr in newCDRRegions:
            if (cdr[2] - cdr[1]) > 5000:
                line = f"{cdr[0]}\t{cdr[1]}\t{cdr[2]}\tCDR\t0\t.\t{cdr[1]}\t{cdr[2]}\t0,25,100\n"
            else:
                line = f"{cdr[0]}\t{cdr[1]}\t{cdr[2]}\tsmall_CDR\t0\t.\t{cdr[1]}\t{cdr[2]}\t0,25,150\n"
            output_lines.append( line )

        for transition in newCDRTransitions:
            line = f"{transition[0]}\t{transition[1]}\t{transition[2]}\tCDR_Transition\t0\t.\t{transition[1]}\t{transition[2]}\t0,25,250\n"
            output_lines.append( line )

        #print( output_lines )

        # Output CDR transition regions to a separate BED file
        with open(outputPrefix, "w") as file:
            for line in output_lines:
                file.write(line)

    def generateLogFiles(self, emissions, states, cpgThresholds, emissionMatrix, transitionMatrix, outputPrefix):
        '''
        Creates log files that contain the emission boundaries + emission/transition matrices

        Attributes:
        - cpgThresholds: 
        - emissionMatrix: 
        - transitionMatrix:
        - outputPrefix:
        '''

        outputPrefix = outputPrefix.replace('.bed', '')
        # write the file with emission boundaries
        thresholds_output = outputPrefix + '.emission_boundaries.csv'
        with open(thresholds_output, "w") as file:
            for threshold in cpgThresholds:
                file.write(str(threshold) + ', ')
        
        emissionMatrix_output = outputPrefix + '.emission_matrix.csv'
        emissionMatrix_df = pd.DataFrame(index=emissions, columns=states)
        for key, value in emissionMatrix.items():
            row, col = key[1], key[0]
            emissionMatrix_df.at[row, col] = value
        emissionMatrix_df.to_csv(emissionMatrix_output)

        transitionMatrix_output = outputPrefix + '.transition_matrix.csv'
        transitionMatrix_df = pd.DataFrame(index=states, columns=states)
        for key, value in transitionMatrix.items():
            row, col = key[1], key[0]
            transitionMatrix_df.at[row, col] = value
        transitionMatrix_df.to_csv(transitionMatrix_output)


def main(options=None):
    '''
    Initializes the first transition and emission matrices and then runs Viterbi Learning to estimate CDR Regions
    '''
    # Creates command line obejct
    thisCommandLine = None
    # Checks if options/parameters were entered
    if options is None:
        thisCommandLine = CommandLine()
    else:
        thisCommandLine = CommandLine(options)
    
    # Creates InitialMatricesEstimate object to get initial transition and emission matrix
    matrixEstimator = InitialMatricesEstimate(thisCommandLine.args.modPosFile, 
                                                thisCommandLine.args.strictCDRFile)
    initialTransitionMatrix = matrixEstimator.getInitialTransitionMatrix()
    path, initialEmissionMatrix = matrixEstimator.getInitialEmissionMatrix()
    chrName = matrixEstimator.getChrName()
    cpgSitesAndProbs = matrixEstimator.getCpGSitesAndProbs()

    # States and emission of HMM Model
    states = ["A", "B", "C"] # A - Normal, B - Transition, C - CDR
    emissions = ["w", "x", "y", "z"] # w - 1/4 quartile of methylation distribution, x - 2/4 quartile of methylation distribution, 
                                        # y - 3/4 quartile of methylation distribution, z - 4/4 quartile of methylation distribution

    # Creates a list for input for the ViterbiLearning object to output optimal CDR regions in a bed file
    viterbiData = [path, emissions, states, 
                   initialTransitionMatrix, initialEmissionMatrix, 
                   float(thisCommandLine.args.learningRate), int(thisCommandLine.args.maxSteps)]
    
    vitLearn = ViterbiLearning(viterbiData)
    estimatedStates = vitLearn.performViterbiLearning()

    print( thisCommandLine.args.outputPrefix )
    print( estimatedStates )

    vitLearn.generateBedFile(chrName, 
                             estimatedStates, 
                             cpgSitesAndProbs, 
                             thisCommandLine.args.outputPrefix)
    
    vitLearn.generateLogFiles(emissions, 
                              states,
                              matrixEstimator.cpgThresholds, 
                              vitLearn.emissionMatrix, 
                              vitLearn.transitionMatrix, 
                              thisCommandLine.args.outputPrefix)


if __name__ == '__main__':
    main()
