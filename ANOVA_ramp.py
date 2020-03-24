import sys
import csv
import numpy as np
import scipy.stats as stats

alphaSub = 0.00015385  # .05/26choose2
alphaSup = 0.005  # .05/5choose2
subpopulations = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB',
                  'GWD', 'ESN', 'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK',
                  'ASW', 'MXL', 'TSI', 'GIH']
superpopulations = ['America', 'East Asia', 'South Asia', 'Africa', 'Europe']


def ANOVAtoCSV(d, csvName):
    header = ['isoform', 'p_value']
    with open(csvName, 'w') as anova_csv_file:
        writer = csv.writer(anova_csv_file)
        writer.writerow(header)
        for iso, p_val in d.items():
            writer.writerow([iso, p_val])


def TukeyToCSV(d, csvName):
    header = ['isoform', 'population1', 'population2', 'p_value']
    with open(csvName, 'w') as tukey_csv_file:
        writer = csv.writer(tukey_csv_file)
        writer.writerow(header)
        for iso, p_val_dict in d.items():
            for p_val, pops in p_val_dict.items():
                pops.sort()
                writer.writerow([iso, pops[0], pops[1], p_val])


def subAnova():
    inputs = []
    for subpopul in subpopulations:
        inputs.append(popDict[subpopul])

    F, p = stats.f_oneway(*inputs)

    return p, tukeyTtest(True)


def superAnova():
    inputs = []
    for superpopul in superpopulations:
        inputs.append(popDict[superpopul])

    F, p = stats.f_oneway(*inputs)

    return p, tukeyTtest(False)


def getAllAnova(sub):
    if sub:
        return subAnova()
    else:
        return superAnova()


def tukeyTtest(sub):
    td = {}
    if sub:
        subpop = [pop for pop in subpopulations]
        for pop1 in subpopulations:
            subpop.pop(0)
            for pop2 in subpop:
                T, P = stats.ttest_ind(np.asarray(popDict[pop1]), np.asarray(popDict[pop2]))
                if P < alphaSub:
                    td[P] = [pop1, pop2]
    else:
        suppop = [pop for pop in superpopulations]
        for pop1 in superpopulations:
            suppop.pop(0)
            for pop2 in suppop:
                T, P = stats.ttest_ind(np.asarray(popDict[pop1]), np.asarray(popDict[pop2]))
                if P < alphaSup:
                    td[P] = [pop1, pop2]
    return td


# aight imma start the program here but maybe pull it out and make more functions? h*ck

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
    popDict = {}
    for superpopulation in superpopulations:
        popDict[superpopulation] = []
    for subpopulation in subpopulations:
        popDict[subpopulation] = []
    with open(sys.argv[1]) as csv_file:  # and arg[1] is the file were going through
        head = csv_file.readline()[:-1].split(",")
        for line in csv_file:
            if line.find(isoformComma) != -1:
                row = line.rstrip().split(",")
                found = True
                popDict[row[1]].append(float(row[6]))
                popDict[row[1]].append(float(row[6]))
    if found:
        subANOVAResults[isoform], subTukeyResults[isoform] = getAllAnova(runForSub)
        superANOVAResults[isoform], superTukeyResults[isoform] = getAllAnova(not runForSub)

ANOVAtoCSV(superANOVAResults, "ramp/" + sys.argv[2] + 'superANOVAResults.csv')
ANOVAtoCSV(subANOVAResults, "ramp/" + sys.argv[2] + 'subANOVAResults.csv')
TukeyToCSV(superTukeyResults, "ramp/" + sys.argv[2] + 'superTukeyResults.csv')
TukeyToCSV(subTukeyResults, "ramp/" + sys.argv[2] + 'subTukeyResults.csv')
