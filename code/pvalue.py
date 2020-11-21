import argparse
import random
from collections import defaultdict
from tanimoto import tanimoto_score as ts


def prot_profile(id_a, id_b, targetfile):
    """
    Given the ID numbers of two proteins, returns two lists of the associated drugs.
    :param id_a: str protein ID number
    :param id_b: str protein ID number
    :param targetfile: read file
    :return: lists
    """
    drugs_a = []
    drugs_b = []
    for line in targetfile:
        line = line.split(',')
        # If the ID matches the protein ID, add the drug to the
        # list of drugs that affect this protein.
        if line[1] == id_a:
            drugs_a.append(line[0])
        # If the ID matches the protein ID, add the drug to the
        # list of drugs that affect this protein.
        elif line[1] == id_b:
            drugs_b.append(line[0])
    return drugs_a, drugs_b


def fingerprint(drugfile):
    """
    Create a dictionary of drugs (str) and their associated fingerprint (set)
    :param drugfile: read file
    :return: dict{str:set}
    """
    dict = defaultdict(set)
    for line in drugfile:
        line = line.split(',')
        dict[line[0]] = set(line[2].split())
    return dict


def t_summary(drugs_a, drugs_b, drug_fpnts):
    """
    Calculates the Tanimoto summary by examining the Tanimoto scores of every drug:drug fingerprint
    combination between two proteins. Sums up only the scores greater than a 0.5 threshold.
    :param drugs_a: list
    :param drugs_b: list
    :param drug_fpnts: dict{str:set}
    :return: float
    """
    # All possible pairs of drugs between the two proteins' profiles
    drug_pairs = [(a, b) for a in drugs_a for b in drugs_b]

    # Calculate T_summary for all Tanimoto sums between two sets of fingerprints
    tanimoto_sum = 0
    for pair in drug_pairs:
        # Get the fingerprint associated with each drug
        da_fpt = drug_fpnts[pair[0]]
        db_fpt = drug_fpnts[pair[1]]
        # Calculate Tanimoto sum
        tanimoto = ts(da_fpt, db_fpt)
        # Include the sum in the summary only if it exceeds 0.5
        if tanimoto > 0.5:
            tanimoto_sum += tanimoto
    return tanimoto_sum


def bootstrap_p(drug_fpnts, drugs_a, drugs_b, r, n, t_sum):
    """
    Find the Tanimoto summary of n randomly chosen drug fingerprints.
    Count how many times these random T summaries exceed the originally calculated T summary.
    Divide by the number of iterations n to ge the bootstrapped p-value.
    :param drug_fpnts: dict{str:set}
    :param drugs_a: list
    :param drugs_b: list
    :param r: int
    :param n: int
    :param t_sum: float
    :return: float
    """
    random.seed(a=r)

    # For n iterations,
    sig_count = 0
    for _ in range(n):
        # Randomly sample w/ replacement n times, where n = number of drugs in original profile
        a_samples = random.choices(list(drug_fpnts.keys()), k=len(drugs_a))
        b_samples = random.choices(list(drug_fpnts.keys()), k=len(drugs_b))

        # Tanimoto score summary for random samples
        sample_summary = t_summary(a_samples, b_samples, drug_fpnts)
        # Count the times the random summary exceeded the measured T summary
        if sample_summary > t_sum:
            sig_count += 1
    return sig_count / n


def main():
    """
    Calculates the bootstrap p-value of a Tanimoto summary score between
    two input proteins.
    :return: int p-value
    """

    parser = argparse.ArgumentParser(description='Generate bootstrap p-value for the comparison of two proteins')

    parser.add_argument('-n', nargs='?', default=500, type=int, help='number of iterations')
    parser.add_argument('-r', nargs='?', default=214, type=int, help='pseudo-random seed value')
    parser.add_argument('drugfile', type=argparse.FileType('r'), help='list of drugs and their fingerprints')
    parser.add_argument('targetfile', type=argparse.FileType('r'), help='drugs and their protein targets')
    parser.add_argument('prot_a', type=str, help='uniprot accession ID')
    parser.add_argument('prot_b', type=str, help='uniprot accession ID')

    args = parser.parse_args()

    # Create a list of drugs that bind each protein
    drugs_a, drugs_b = prot_profile(args.prot_a, args.prot_b, args.targetfile)

    # Dictionary: {drug: fingerprint}
    drug_fpnts = fingerprint(args.drugfile)

    # Calculate the Tanimoto summary between the two drug profiles
    t_sum = t_summary(drugs_a, drugs_b, drug_fpnts)

    # Determine the p-value of t_sum by random sampling
    p_value = bootstrap_p(drug_fpnts, drugs_a, drugs_b, args.r, args.n, t_sum)

    print(p_value)
    return p_value


if __name__ == '__main__':
    main()
