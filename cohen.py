import sys
import csv
from numpy import std, mean, sqrt
import math

codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
          'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
          'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
          'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT',
          'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
subpopulations = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB',
                  'GWD', 'ESN', 'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK',
                  'ASW', 'MXL', 'TSI', 'GIH']
superpopulations = ['America', 'East Asia', 'South Asia', 'Africa', 'Europe']


def cohenToCSV(d, csvName):
    header = ['isoform', 'codon', 'population1', 'population2', 'cohens_d']
    with open(csvName, 'w') as cohen_csv_file:
        writer = csv.writer(cohen_csv_file)
        writer.writerow(header)
        for iso, codonDict in d.items():
            for c, tup_list in codonDict.items():
                for tup in tup_list:
                    writer.writerow([iso, c, tup[0], tup[1], str(tup[2])])


def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    D = (mean(x) - mean(y)) / sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)
    if math.isinf(D):
        D = "D test broke"
    return D


def cohenTest(r):
    pop1List = isoDict[r[0]][r[1]][r[2]]
    pop2List = isoDict[r[0]][r[1]][r[3]]
    return cohen_d(pop1List, pop2List)


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
with open(sys.argv[2]) as csv_file:
    for line in csv_file:
        row = line.rstrip().split(",")
        isoforms.append(row[1])

# initialize results dicts

superDResults = {}
subDResults = {}
isoDict = {}
runForSub = True

# fill dictionaries
for isoform in isoforms:
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
        isoDict[isoform] = codonDicts

with open(sys.argv[3]) as csv_file:
    csv_file.readline()
    for line in csv_file:
        Row = line.rstrip().split(",")
        if Row[0] in isoDict.keys():
            if Row[0] in superDResults.keys():
                if Row[1] in superDResults[Row[0]].keys():
                    superDResults[Row[0]][Row[1]].append((Row[2], Row[3], cohenTest(Row)))
                else:
                    superDResults[Row[0]][Row[1]] = [(Row[2], Row[3], cohenTest(Row))]
            else:
                superDResults[Row[0]] = {}
                superDResults[Row[0]][Row[1]] = [(Row[2], Row[3], cohenTest(Row))]

with open(sys.argv[4]) as csv_file:
    csv_file.readline()
    for line in csv_file:
        Row = line.rstrip().split(",")
        if Row[0] in isoDict.keys():
            if Row[0] in subDResults.keys():
                if Row[1] in subDResults[Row[0]].keys():
                    subDResults[Row[0]][Row[1]].append((Row[2], Row[3], cohenTest(Row)))
                else:
                    subDResults[Row[0]][Row[1]] = []
                    subDResults[Row[0]][Row[1]].append((Row[2], Row[3], cohenTest(Row)))
            else:
                subDResults[Row[0]] = {}
                subDResults[Row[0]][Row[1]] = []
                subDResults[Row[0]][Row[1]].append((Row[2], Row[3], cohenTest(Row)))

cohenToCSV(superDResults, 'codon_counts/' + sys.argv[2] + 'superDResults.csv')
cohenToCSV(subDResults, 'codon_counts/' + sys.argv[2] + 'subDResults.csv')
