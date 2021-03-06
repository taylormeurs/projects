import sys
import csv
import numpy as np
import scipy.stats as stats

alphaSub = 0.00015385  # .05/26choose2
alphaSup = 0.005  # .05/5choose2
codons = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Z']
subpopulations = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB',
                  'GWD', 'ESN', 'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK',
                  'ASW', 'MXL', 'TSI', 'GIH']
superpopulations = ['America', 'East Asia', 'South Asia', 'Africa', 'Europe']


def ANOVAtoCSV(d, csvName):
    header = ['isoform', 'amino acid', 'p_value']
    with open(csvName, 'w') as anova_csv_file:
        writer = csv.writer(anova_csv_file)
        writer.writerow(header)
        for iso, codonDict in d.items():
            for c, p_val in codonDict.items():
                writer.writerow([iso, c, p_val])


def TukeyToCSV(d, csvName):
    header = ['isoform', 'amino acid', 'population1', 'population2', 'p_value']
    with open(csvName, 'w') as tukey_csv_file:
        writer = csv.writer(tukey_csv_file)
        writer.writerow(header)
        for iso, codonDict in d.items():
            for c, p_val_dict in codonDict.items():
                for p_val, pops in p_val_dict.items():
                    pops.sort()
                    writer.writerow([iso, c, pops[0], pops[1], p_val])


def subAnova():
    resultsDict = {}
    tukeyDict = {}

    for curCodon in codons:
        inputs = []
        for subpopul in subpopulations:
            inputs.append(codonDicts[curCodon][subpopul])

        F, p = stats.f_oneway(*inputs)
        if p < alphaSub:
            resultsDict[curCodon] = p
            tukeyDict[curCodon] = tukeyTtest(codonDicts[curCodon], True)

    return resultsDict, tukeyDict


def superAnova():
    resultsDict = {}
    tukeyDict = {}

    for curCodon in codons:
        inputs = []
        for superpopul in superpopulations:
            inputs.append(codonDicts[curCodon][superpopul])

        F, p = stats.f_oneway(*inputs)
        if p < alphaSup:
            resultsDict[curCodon] = p
            tukeyDict[curCodon] = tukeyTtest(codonDicts[curCodon], False)

    return resultsDict, tukeyDict


def getAllAnova(sub):
    if sub:
        return subAnova()
    else:
        return superAnova()


def tukeyTtest(curDict, sub):
    td = {}
    if sub:
        subpop = [pop for pop in subpopulations]
        for pop1 in subpopulations:
            subpop.pop(0)
            for pop2 in subpop:
                T, P = stats.ttest_ind(np.asarray(curDict[pop1]), np.asarray(curDict[pop2]))
                if P < alphaSub:
                    td[P] = [pop1, pop2]        #todo if pvals are exactly the same it deletes one...
    else:
        suppop = [pop for pop in superpopulations]
        for pop1 in superpopulations:
            suppop.pop(0)
            for pop2 in suppop:
                T, P = stats.ttest_ind(np.asarray(curDict[pop1]), np.asarray(curDict[pop2]))
                if P < alphaSup:
                    td[P] = [pop1, pop2]
    return td


# aight imma start the program here but maybe pull it out and make more functions? h*ck

# make population map
subpopMap = {}
superpopMap = {}
with open('sampleKeyPop.csv') as csv_file:
    for line in csv_file:
        row = line.rstrip().split(",")
        subpopMap[row[0]] = row[1]
        superpopMap[row[0]] = row[2]

# get list of isoforms from sys.argv[2]
isoforms = []
with open("splitIsoforms/" + sys.argv[2]) as csv_file:
    for line in csv_file:
        row = line.rstrip().split(",")
        isoforms.append(row[1])

# initialize results dicts

subANOVAResults = {}
subTukeyResults = {}
superANOVAResults = {}
superTukeyResults = {}
runForSub = True

# fill dictionaries
for isoform in isoforms:
    print(isoform)
    isoformComma = "," + isoform + ","
    found = False
    # reinitialize codonDicts
    codonDicts = {}
    for codon in codons:
        codonDicts[codon] = {}
        for superpopulation in superpopulations:
            codonDicts[codon][superpopulation] = []
        for subpopulation in subpopulations:
            codonDicts[codon][subpopulation] = []
    with open(sys.argv[1]) as csv_file:  # and arg[1] is the file were going through
        head = csv_file.readline()[:-1].split(",")
        for line in csv_file:
            if line.find(isoformComma) != -1:
                row = line.rstrip().split(",")
                found = True
                for codon in codons:
                    codonDicts[codon][subpopMap[row[0]]].append(int(row[head.index(codon)]))
                    codonDicts[codon][superpopMap[row[0]]].append(int(row[head.index(codon)]))
    if found:
        subANOVAResults[isoform], subTukeyResults[isoform] = getAllAnova(runForSub)
        superANOVAResults[isoform], superTukeyResults[isoform] = getAllAnova(not runForSub)

ANOVAtoCSV(superANOVAResults, "co_pair/" + sys.argv[2] + 'superANOVAResults.csv')
ANOVAtoCSV(subANOVAResults, "co_pair/" + sys.argv[2] + 'subANOVAResults.csv')
TukeyToCSV(superTukeyResults, "co_pair/" + sys.argv[2] + 'superTukeyResults.csv')
TukeyToCSV(subTukeyResults, "co_pair/" + sys.argv[2] + 'subTukeyResults.csv')
