import sys
import itertools
import pvalue


def main():
    """
    Creates a file with the pairs of proteins which have demonstrated
    a p-value of significance, which can be used as nodes in network generation.
    :return: txt file
    """
    # For all the combinations of protein pairs within the protein_nodes file,
        # If the bootstrap p-value is <= 0.05,
            # Write the pairs' uniprot accension ID's to the output file.

    if len(sys.argv) < 4:
        print('Please specify drug file, target file, and protein node file.')

    drugfile = sys.argv[1]
    targetfile = sys.argv[2]
    protein_node = sys.argv[3]

    # Default values for random seed (r) and iterations for bootstrapping (n)
    r = 214
    n = 500
    sig_p = 0.05

    sig_prots = set()

    with open(protein_node, 'r') as f:
        prots = f.readlines()[1:]

    with open(drugfile, 'r') as d:
        drug_fpnts = pvalue.fingerprint(d)

    with open('output/network_edgelist.txt', 'w') as e:
        for a, b in itertools.combinations(prots, r=2):
            # Parse the protein ID numbers
            accession_a = a.split(',')[0]
            accession_b = b.split(',')[0]

            with open(targetfile, 'r') as t:
                # Use the protein ID numbers to gather the list of associated drugs
                drugs_a, drugs_b = pvalue.prot_profile(accession_a, accession_b, t)
            # Calculate Tanimoto summary for the proteins by looking at
            # Tanimoto score of drug fingerprint pairs
            t_sum = pvalue.t_summary(drugs_a, drugs_b, drug_fpnts)

            if pvalue.bootstrap_p(drug_fpnts, drugs_a, drugs_b, r, n, t_sum) <= sig_p:
                e.write(f'{accession_a} {accession_b}\n')

                sig_prots.add(accession_a)
                sig_prots.add(accession_b)


if __name__ == '__main__':
    main()
