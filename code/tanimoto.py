import itertools
import sys
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import math


def parse_target(targetfile):
    """
    Read target file, create a dict with drugs and their target proteins.
    :param targetfile: file path
    :return: dict {drug: protein}
    """
    drug_prot_dict = defaultdict(list)
    with open(targetfile, 'r') as j:
        # Eliminate header
        targets = j.readlines()[1:]
        # Create a list of protein targets for each drug
        for row in targets:
            row = row.split(sep=',')
            drug_id = row[0]
            prot_id = row[1]
            drug_prot_dict[drug_id].append(prot_id)
    return drug_prot_dict


def histo(scores, type):
    """
    Saves a histogram of Tanimoto summary scores from 1 of 3 groups:
    1. Shares a target protein
    2. Does not share a target protein
    3. All scores regardless of protein
    :param scores: list of drugs belonging to respective group
    :param type: str group description
    :return: saves plot
    """
    # My personal SUNet
    sunet = 'seshwan2'
    # Two different bin parameter strategies to use: rice or square root
    sqr_bins = round(math.sqrt(len(scores)))
    bins_rice = round(2*(len(scores)**(1./3)))

    plt.figure()
    # Rice method looks a bit better for this data
    plt.hist(scores, bins=bins_rice)
    plt.title(sunet + ' ' + type)
    plt.ylabel('Frequency')
    plt.xlabel('Tanimoto Score')
    plt.savefig('output/{}_tanimoto.png'.format(type))


def fingerprint_id(line):
    """
    Parse out the drug ID number and its corresponding fingerprint
    :param line: list of drug ID, generic name, and fingerprint
    :return: str id
    :return: set fingerprint
    """
    line = line.split(',')
    id = line[0]
    fingerprint = set(line[2].split())
    return id, fingerprint


def tanimoto_score(fpnt1, fpnt2):
    """
    Calculates the tanimoto score between two sets of drug fingerprints
    :param fpnt1: set of fingerprints
    :param fpnt2: set of fingerprints
    :return: float
    """
    # Number of elements in union
    shared = len(fpnt1 & fpnt2)
    # Number of unique elements between both
    union = len(fpnt1 | fpnt2)
    # Tanimoto
    return float(shared / union)


def main():
    """
    Find the Tanimoto coefficient between drugs based on fingerprint data.
    Determine whether or not drugs have any target proteins in common.
    :return: 3 histograms
    """

    if len(sys.argv) != 4:
        print('Please specify drugs file, targets file, and output file')

    drugfile = sys.argv[1]
    targetfile = sys.argv[2]
    outputfile = sys.argv[3]

    drug_prot_dict = parse_target(targetfile)

    # Read drug file
    with open(drugfile, 'r') as f:
        drugs = f.readlines()[1:]

    # Inputs for graph function
    total = []
    shared = []
    not_shared = []

    # To be used to output 6 decimal places in the end
    epsilon = 0.000000001
    # Write output file
    with open(outputfile, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        # For each unique combination of drugs,
        for (first, second) in itertools.combinations(drugs, r=2):
            # Parse out the line to extract the ID and macc
            id_1, fingerprint_1 = fingerprint_id(first)
            id_2, fingerprint_2 = fingerprint_id(second)
            # Run prints through Tanimoto function
            t_score = tanimoto_score(fingerprint_1, fingerprint_2) + epsilon
            total.append(t_score)
            # Check to see if the two drugs share any target proteins
            if any(item in drug_prot_dict[id_1] for item in drug_prot_dict[id_2]):
                common = 1
                shared.append(t_score)
            else:
                common = 0
                not_shared.append(t_score)
            # Write output
            writer.writerow([id_1, id_2, '{:.6f}'.format(t_score), common])

    # Create histogram of 3 different groups of data
    histo(total, 'all')
    histo(shared, 'shared')
    histo(not_shared, 'notshared')


if __name__ == '__main__':
    main()
