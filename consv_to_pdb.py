import math
import matplotlib.pyplot as plt
import numpy as np
from tkinter import filedialog

# AA background frequencies (https://www.ebi.ac.uk/uniprot/TrEMBLstats)
back_freq = {
    'A': 0.0913,
    'R': 0.0580,
    'N': 0.0380,
    'D': 0.0547,
    'C': 0.0127,
    'E': 0.0618,
    'Q': 0.0377,
    'G': 0.0730,
    'H': 0.0220,
    'I': 0.0556,
    'L': 0.0988,
    'K': 0.0499,
    'M': 0.0234,
    'F': 0.0390,
    'P': 0.0492,
    'S': 0.0671,
    'T': 0.0557,
    'W': 0.0130,
    'Y': 0.0290,
    'V': 0.0693,
    'X': 1.0,
    '-': 1.0,
    'B': 0.0001
}

# def readfile(filename):
#     with open(filename) as file:
#
#         seqs = [part[2].replace('\n', '') for part in [entry.partition('\n') for entry in file.read().split('>')[1:]]]
#     return seqs

def readfile(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Extract sequences and their headers
    seqs = []
    headers = []
    missing_map = []
    seqs_dict = {}
    for line in lines:
        if line.startswith('>'):
            headers.append(line.strip())
            seqs.append('')
        else:
            seqs[-1] += line.strip()

    for i in range(len(headers)):
        number_map = {}
        number_map2 = {}
        bar_count = 0
        for j in range(len(seqs[i])):
            if seqs[i][j] == '-':
                bar_count += 1
                missing_map.append(j + 1)
            number_map[j + 1 - bar_count] = j+1  # number_before_msa = key, number_after_msa = value
            number_map2[j+1] = j + 1 - bar_count # number_after_msa = key
        seqs_dict[headers[i].split('|')[1]] = seqs[i], number_map, number_map2  # key = accessions

    print(seqs_dict)
    missing_map = list(set(missing_map))
    print(missing_map)
    return seqs, headers, seqs_dict, missing_map

def getcolumns(seqs):
    n = len(seqs[0])
    columns = []
    AAtype = []

    for j in range(n):
        characters = []
        for i in range(len(seqs)):
            characters.append(seqs[i][j])
        columns.append(characters)
        AAtype.append(list(set(columns[j])))

    return columns, AAtype

def freqs(AAtype, columns, back_freq):
    freqList = []
    adjFList = []

    for i in range(len(columns)):
        freq = []
        adjF = []
        for symbol in AAtype[i]:
            ctr = columns[i].count(symbol)
            charFreq = float(ctr) / len(columns[i])
            freq.append(charFreq)
            adjF.append(charFreq / back_freq[symbol])
        freqList.append(freq)
        adjFList.append(adjF)

    return AAtype, freqList, adjFList

def shannon(AAtype, freqList, adjFList,columns):
    shan_ent = []
    shan_ent2 = []

    for i in range(len(columns)):
        charEnt = []
        charEnt2 = []

        for j in range(len(freqList[i])):
            charEnt.append(freqList[i][j] * math.log2(adjFList[i][j]))
            charEnt2.append((freqList[i][j] * math.log2(adjFList[i][j])) ** 2)

        shan_ent.append(-sum(charEnt))
        shan_ent2.append(sum(charEnt2) / len(AAtype[i]))

    return shan_ent2

def plot_conserv(shan_ent2):
    conserv = []
    for i in range(len(shan_ent2)):
        # print(i + 1, shan_ent2[i] ** 0.5)
        conserv.append(shan_ent2[i] ** 0.5)

    plt.plot(range(len(conserv)), conserv)
    plt.title('Position-specific conservation')
    plt.xlabel('Position')
    plt.ylabel('Conservation')
    # plt.show()

def get_conservation(seqs):
    columns, AAtype = getcolumns(seqs)
    AAtype, freqList, adjFList = freqs(AAtype, columns, back_freq)
    entropy = shannon(AAtype, freqList, adjFList, columns)
    plot_conserv(entropy)


def main():
    filename = filedialog.askopenfilename(title='Select the file you want to use as the reference')
    seqs, headers, seqs_dictionary, missing_map = readfile(filename)
    get_conservation(seqs)


# main()
