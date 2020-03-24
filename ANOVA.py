import sys
import csv
import numpy as np
import scipy.stats as stats

alphaSub = 0.00015385  # .05/26choose2
alphaSup = 0.005  # .05/5choose2
codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
          'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
          'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
          'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT',
          'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
subpopulations = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB',
                  'GWD', 'ESN', 'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK',
                  'ASW', 'MXL', 'TSI', 'GIH']
superpopulations = ['America', 'East Asia', 'South Asia', 'Africa', 'Europe']


def ANOVAtoCSV(d, csvName):
    header = ['isoform', 'codon', 'p_value']
    with open(csvName, 'w') as anova_csv_file:
        writer = csv.writer(anova_csv_file)
        writer.writerow(header)
        for iso, codonDict in d.items():
            for c, p_val in codonDict.items():
                writer.writerow([iso, c, p_val])


def TukeyToCSV(d, csvName):
    header = ['isoform', 'codon', 'population1', 'population2', 'p_value']
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
                print(isoform, P)
                if P < alphaSub:
                    td[P] = [pop1, pop2]
    else:
        suppop = [pop for pop in superpopulations]
        for pop1 in superpopulations:
            suppop.pop(0)
            for pop2 in suppop:
                T, P = stats.ttest_ind(np.asarray(curDict[pop1]), np.asarray(curDict[pop2]))
                print(isoform, P)
                if P < alphaSup:
                    td[P] = [pop1, pop2]
    return td


# aight imma start the program here but maybe pull it out and make more functions? h*ck

# make population map
subpopMap = {}
superpopMap = {}
with open('sampleKeyPop.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        subpopMap[row[0]] = row[1]
        superpopMap[row[0]] = row[2]

# get list of isoforms
isoforms = []
with open('longestIsoforms') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        isoforms.append(row[1])

# initialize results dicts

subANOVAResults = {}
subTukeyResults = {}
superANOVAResults = {}
superTukeyResults = {}
runForSub = True

# fill dictionaries
for isoform in isoforms:
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
        csv_reader = csv.reader(csv_file, delimiter=',')
        head = next(csv_reader)
        for row in csv_reader:
            if row[2] == isoform:
                found = True
                for codon in codons:
                    codonDicts[codon][subpopMap[row[0]]].append(int(row[head.index(codon)]))
                    codonDicts[codon][superpopMap[row[0]]].append(int(row[head.index(codon)]))
    if found:
        subANOVAResults[isoform], subTukeyResults[isoform] = getAllAnova(runForSub)
        superANOVAResults[isoform], superTukeyResults[isoform] = getAllAnova(not runForSub)

ANOVAtoCSV(superANOVAResults, 'superANOVAResults.csv')
ANOVAtoCSV(subANOVAResults, 'subANOVAResults.csv')
TukeyToCSV(superTukeyResults, 'superTukeyResults.csv')
TukeyToCSV(subTukeyResults, 'subTukeyResults.csv')
