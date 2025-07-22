def parse_interactome(interactome_file):
    '''
    Parses an interactome SIF file.

    Arguments:
    - interactome_file: str, path to SIF file (3 tab-separated columns: ENSG1 pp ENSG2)

    Returns:
    - nbNodes: int
    - weights: list of floats, flattened adjacency matrix,
      weights[i*nbCols + j] contains the weight from node i to node j,
      note: weights must be 0/1 for unweighted, or in [0, 1] for weighted
    '''
    nodes = set()
    edges = []

    try:
        with open(interactome_file, 'r') as f:
            for line in f:
                split_line = line.rstrip().split('\t')
                if len(split_line) != 3:
                    raise Exception("Bad line in the interactome file")

                gene1, weight, gene2 = split_line
                if weight != "pp":
                    try:
                        weight = float(weight)
                    except ValueError:
                        raise Exception(f"Invalid weight '{weight}' on line: {line}")
                else:
                    weight = 1.0

                edges.append((gene1, gene2, weight))
                nodes.update([gene1, gene2])

    except Exception as e:
        raise Exception(f"Cannot open provided interactome file {e}")

    nodes = sorted(nodes)
    nbNodes = len(nodes)
    node_index = {gene: i for i, gene in enumerate(nodes)}

    # initialize weights
    weights = [0.0] * (nbNodes * nbNodes)

    for gene1, gene2, weight in edges:
        i, j = node_index[gene1], node_index[gene2]
        weights[i * nbNodes + j] = weight

    return nbNodes, weights
