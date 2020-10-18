


def read_fastq(fichier : str):
    with open(fichier) as filin:
        for line in enumerate(filin):
            yield next(filin)[:-1]
            next(filin)
            next(filin)

def cut_kmer(seq : str,k : int):
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

def build_kmer_dict(fichier : str, k : int):
    out = {}
    for i in read_fastq(fichier):
        for j in cut_kmer(i,k = k ):
            if not j in out :
                out.setdefault(j, 1)
            else :
                out[j] += 1
    return out



if __name__ == "__main__" :
    k_dict = build_kmer_dict("../data/eva71_plus_perfect.fq", k = 10)
    


