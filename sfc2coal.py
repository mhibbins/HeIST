import sys
import re
import numpy as np
import numpy.polynomial.polynomial as poly

def read_newick(file): 
        
        with open(sys.argv[1]) as treefile:
                newick_string = treefile.readline()
        
        return newick_string

def subs2coal(newick_string):
        
        '''
        Takes a newick string with the nodes labelled with concordance factors,
        and returns the same string with branch lengths converted to coalescent
        units.
        '''

        scfs = re.findall("\)(.*?)\:", newick_string) #regex to get concordance factors
        scfs = list(filter(None, scfs)) #removes empty values (root has no scf)
        coal_internals = []

        for i in range(len(scfs)):

                if (float(scfs[i])/100) < 1 and float(scfs[i]) != 1: #catch for missing data / values of 1
                        coal_internals.append(-1*(np.log(1-(float(scfs[i])/100))))
                else:
                        coal_internals.append(np.NaN)

        coal_internals = [float(i) for i in coal_internals]
        newick_internals = []
        internal_sections = []
        sections = newick_string.split(':')

        for i in range(len(sections)):
                if len(sections[i].split(")")) > 1 and sections[i].split(")")[1] in scfs:
                                newick_internals.append(newick_string.split(':')[i+1].split(',')[0].split(")")[0])
                                internal_sections.append(sections[i+1])

        
        newick_internals = [float(i) for i in newick_internals]
         
        def branch_regression(newick_branches, coal_internals):

                '''
                Returns the regression coefficient of the coalescent branch lengths
                on the branch lengths in the newick string
                '''

                for i in range(len(newick_branches)): #drops missing data
                        if np.isnan(newick_branches[i]) or np.isnan(coal_internals[i]):
                                del newick_branches[i]
                                del coal_internals[i]

                intercept, slope = poly.polyfit(newick_branches, coal_internals, 1)

                return intercept, slope

        intercept, coef = branch_regression(newick_internals, coal_internals)
        
        tip_sections = []

        for i in range(len(sections)): #gets the sections with tip lengths
                if sections[i] not in internal_sections:
                        tip_sections.append(sections[i])

        newick_tips = []

        for i in range(len(tip_sections)): #gets the tip lengths
                if len(tip_sections[i].split(',')) > 1:
                        newick_tips.append(tip_sections[i].split(',')[0])
                elif len(tip_sections[i].split(')')) > 1:
                        newick_tips.append(tip_sections[i].split(')')[0])
        
        newick_tips = [float(i) for i in newick_tips]

        coal_tips = [(newick_tips[i]*coef + intercept) for i in range(len(newick_tips))]

        print(len(newick_tips))
        print(len(coal_tips))

        coal_internals = [(newick_internals[i]*coef + intercept) if np.isnan(coal_internals[i]) else coal_internals[i] for i in range(len(newick_internals))]

        lengths = re.findall("\d+\.\d+", newick_string)

        lengths = [lengths[i] for i in range(len(lengths)) if float(lengths[i]) < 1]

        coal_lengths = []

        for i in range(len(lengths)):
                if newick_internals.count(lengths[i]) > 1 or newick_tips.count(lengths[i]) > 1:
                        sys.exit('Error: Duplicate branch lengths')
                elif float(lengths[i]) in newick_internals:
                        coal_lengths.append(coal_internals[newick_internals.index(float(lengths[i]))])
                elif float(lengths[i]) in newick_tips:
                        coal_lengths.append(coal_tips[newick_tips.index(float(lengths[i]))])
                elif float(lengths[i]) == 0: #deals with roots of length 0 in smoothed trees
                        coal_lengths.append(float(0))
        
        coal_lengths = [str(coal_lengths[i]) for i in range(len(coal_lengths))]
       
        coal_newick_string = newick_string

        for i in range(len(lengths)):
             coal_newick_string = coal_newick_string.replace(str(lengths[i]), str(coal_lengths[i]))
    
        for i in range(len(scfs)):
                coal_newick_string = coal_newick_string.replace(str(scfs[i]), '')    
        
        return coal_newick_string

newick_string = read_newick(sys.argv[1])
coal_newick_string = subs2coal(newick_string)
print(newick_string)
print(coal_newick_string)


