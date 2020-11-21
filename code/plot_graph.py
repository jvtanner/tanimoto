import networkx as nx
import matplotlib.pyplot as plt
import sys


def labels_colors(prot_nodes):
    """
    Using prot_nodes file, create a dict linking the accession number to the uniprot ID,
    and another dict linking the uniprot ID to color, based on indication.
    :param prot_nodes: csv file
    :return: dict, dict
    """
    # Maps colors to indications listed in prot_nodes.csv file
    color_codes = {'bp': 'red', 'bp;cholesterol': 'green', 'bp;cholesterol;diabetes': 'blue', 'bp;diabetes': 'purple'}

    relabel = {}
    recolor = {}

    f = open(prot_nodes, 'r')
    header = True
    for line in f.readlines():
        line = line.split(',')
        # Don't read the header
        if header:
            header = False
            continue
        indication = line[2].strip()
        # Link color to indication
        color = color_codes[indication]
        # Create dict {Uniprot_ID: color}
        recolor[line[1]] = color
        # Create relabel dict {Accession_ID: Uniprot_ID}
        relabel[line[0]] = line[1]
    f.close()

    return relabel, recolor


def main():
    """
    Create a visual network of the edgemap created in plot_graph.py.
    Proteins are relabeled with their uniprot ID's.
    Node colors correlate to indication.
    :return: png file of network map
    """
    edgelist = sys.argv[1]
    prot_nodes = sys.argv[2]

    # Two dicts with appropriate relabeling and recoloring schema
    relabel, recolor = labels_colors(prot_nodes)

    # Create plot
    plt.figure(figsize=(8, 8), dpi=150)
    plt.title("Proteins with Significant Tanimoto Summaries")

    # Pass edge list into networkx object
    edge_list = open(edgelist, "rb")
    G = nx.read_edgelist(edge_list)
    edge_list.close()
    # Relabel nodes to uniprot ID using dict created above.
    H = nx.relabel_nodes(G, relabel)
    # Use color dict to create list of color which mirrors the order of nodes.
    color_list = []
    for node in H.nodes():
        color_list.append(recolor[node])
    # Draw network with appropriate colors
    nx.draw(H, node_color=color_list, with_labels=1)
    plt.savefig("output/network.png", format="PNG")
    plt.show()


if __name__ == '__main__':
    main()
