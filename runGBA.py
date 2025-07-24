import ctypes


def GBA_centrality():

    class adjacencyMatrix(ctypes.Structure):
        _fields_ = [('nbCols', ctypes.c_uint),
                    ('weights', ctypes.POINTER(ctypes.c_float))]

    class geneScores(ctypes.Structure):
        _fields_ = [('nbGenes', ctypes.c_uint),
                    ('scores', ctypes.POINTER(ctypes.c_float))]

    so_file = "./gbaCentrality.so"
    gbaLibrary = ctypes.CDLL(so_file)
    # declare function signature
    gbaLibrary.gbaCentrality.argtypes = [
        ctypes.POINTER(adjacencyMatrix),
        ctypes.POINTER(geneScores),
        ctypes.c_float,
        ctypes.POINTER(geneScores)
    ]
    gbaLibrary.gbaCentrality.restype = None

    # generate structures: adjacencyMatrix, geneScores
    weights = [0.00, 0.20, 0.50, 0.00, 0.40, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 1.00, 0.00, 1.00, 1.00, 1.00]
    weights_ctype = (ctypes.c_float * len(weights))(*weights)
    A = adjacencyMatrix(nbCols=4,
                        weights=weights_ctype)

    causal_genes = [1.0, 0.0, 0.0, 1.0]
    causal_genes_ctype = (ctypes.c_float * len(causal_genes))(*causal_genes)
    causal = geneScores(nbGenes=4,
                        scores=causal_genes_ctype)

    scores = (ctypes.c_float * 4)()
    res = geneScores(nbGenes=4, scores=ctypes.cast(scores, ctypes.POINTER(ctypes.c_float)))

    gbaLibrary.gbaCentrality(
        ctypes.byref(A),
        ctypes.byref(causal),
        ctypes.c_float(0.5),
        ctypes.byref(res)
    )

    print("Output scores:", [res.scores[i] for i in range(res.nbGenes)])


if __name__ == "__main__":
    GBA_centrality()
