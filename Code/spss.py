fasta_file = "data/reads_05x.fasta"
complement_cache = {}

def rev_complement(kmer):
    """
    Return the reverse complement of a DNA k-mer.

    This function takes a string representing a k-mer (a nucleotide sequence)
    and returns its reverse complement. The cache is used to optimize repeated
    calls with the same k-mers by avoiding recomputation.

    Parameters:
    kmer (str): The string representing the k-mer, composed of 'A', 'T', 'C', 'G'.

    Returns:
    str: The reverse complement of the k-mer.

    Example:
    >>> rev_complement("ATCG")
    "CGAT"

    Note:
    - The `complement_cache` cache stores previously computed results to improve performance.

    """
    if kmer in complement_cache:
        return complement_cache[kmer]
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rev_comp = ''.join(complement.get(c, '') for c in reversed(kmer))
    complement_cache[kmer] = rev_comp
    return rev_comp

def rev_complement(kmer):
    """
    Return the reverse complement of a DNA k-mer.

    This function takes a string representing a k-mer (a nucleotide sequence)
    and returns its reverse complement.

    Parameters:
    kmer (str): The string representing the k-mer, composed of 'A', 'T', 'C', 'G'.

    Returns:
    str: The reverse complement of the k-mer.

    Example:
    >>> rev_complement("ATCG")
    "CGAT"

    """
    complement = str.maketrans("ATCG", "TAGC")
    return kmer.translate(complement)[::-1]

def canonique(kmer):
    """
    Compute the canonical form of a k-mer.

    The canonical form is the lexicographically smallest string between
    the k-mer and its reverse complement.

    Args:
        kmer (str): The input k-mer.

    Returns:
        str: The canonical k-mer.
    """
    return min(kmer, rev_complement(kmer))

def listing_kmer(fasta_file,k):
    """
    Extract and count k-mers from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.
        k (int): Length of k-mers.

    Returns:
        dict: A dictionary where keys are canonical k-mers and values
              are their abundances (frequencies).
    """
    kmer_dict = {}
    with open(fasta_file,"r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                for num in range(len(line)-k + 1):
                    kmer = canonique(line[num:k + num])
                    if kmer in kmer_dict:
                        kmer_dict[kmer] += 1
                    else :
                        kmer_dict[kmer] = 1
    return kmer_dict

def solide(kmer_dict,t):
    """
    Filter solid k-mers.

    A k-mer is considered "solid" if its abundance is above
    a given threshold.

    Args:
        kmer_dict (dict): Dictionary of k-mers with their abundances.
        t (int): Minimum abundance threshold.

    Returns:
        dict: A dictionary containing only solid k-mers.
    """
    for key in list(kmer_dict.keys()):
        if kmer_dict[key] <= t:
            del kmer_dict[key]
    return kmer_dict

def kmer_histo(dbg, ab_max = 100):
    """
    Build a histogram of k-mer abundances.

    Args:
        dbg (dict): A dictionary containing k-mers and their abundances.
        ab_max (int): Maximum abundance threshold to consider.

    Returns:
        list: A list where each index represents an abundance, and the value
              represents the number of k-mers with that abundance.
    """
    histo = [0 for i in range(ab_max)]
    for kmer in dbg:
        abundance = dbg[kmer]
        if abundance < ab_max:
            histo[abundance] += 1
        else:
            histo[ab_max-1] += 1
    return histo

def right_neighbors(kmer, solide_kmer_dict):
    """
    Find right neighbors of a given k-mer.

    Args:
        kmer (str): The input k-mer.
        solide_kmer_dict (dict): Dictionary of solid k-mers.

    Returns:
        list: A list of valid right-neighbor k-mers.
    """
    res = []
    for nuc in "ATCG":
        voisin = kmer[1:] + nuc
        voisin_canonique = canonique(voisin)
        if voisin_canonique in solide_kmer_dict and solide_kmer_dict[voisin_canonique] > 0:
            res.append(voisin)
    return res

def right_extend( kmer, solide_kmer_dict):
    """
    Extend a k-mer to the right by adding a nucleotide.

    Args:
        kmer (str): The k-mer to extend.
        solide_kmer_dict (dict): Dictionary of solid k-mers.

    Returns:
        str: The added nucleotide.
    """
    nuc = right_neighbors(kmer,solide_kmer_dict)[0][-1]
    solide_kmer_dict[kmer] = 0
    return nuc

def check_rev_comp(kmer,solide_kmer_dict):
    """
    Check whether the reverse complement of a k-mer is solid.

    Args:
        kmer (str): The input k-mer.

    Returns:
        bool: True if the reverse complement is solid, False otherwise.
    """
    kmer_canonique = canonique(kmer)
    return kmer_canonique in solide_kmer_dict and solide_kmer_dict[kmer_canonique] > 0

def extend(solide_kmer_dict):
    """
    Extend solid k-mers to build the SPSS.

    Args:
        solide_kmer_dict (dict): Dictionary of solid k-mers.

    Returns:
        list: The SPSS.
    """
    res = []
    for kmer in list(solide_kmer_dict.keys()): 
        voisin = kmer
        if solide_kmer_dict[voisin] != 0:
            kmer_actu = ""
            chemin = kmer
            while check_rev_comp(voisin,solide_kmer_dict) and len(right_neighbors(voisin,solide_kmer_dict)) > 0:
                nuc_ajoute = right_extend(voisin,solide_kmer_dict)
                kmer_actu = voisin
                voisin = kmer_actu[1:] + nuc_ajoute
                solide_kmer_dict[canonique(kmer_actu)] = 0
                chemin += nuc_ajoute
            kmer_actu = ""
            while check_rev_comp(voisin,solide_kmer_dict) and len(right_neighbors(rev_complement(voisin),solide_kmer_dict)) > 0:
                nuc_ajoute = rev_complement(right_extend(rev_complement(voisin),solide_kmer_dict))
                kmer_actu = voisin
                voisin = nuc_ajoute + voisin[:-1]
                solide_kmer_dict[canonique(kmer_actu)] = 0
                chemin = nuc_ajoute + chemin
            res.append(chemin)
            solide_kmer_dict[canonique(voisin)] = 0
    return res

def concat_spss(spss):
    """
    Concatenate reconstructed sequences into a single string.

    Sequences are separated by the '%' character and the string ends with '$'.

    Args:
        spss (list): List of reconstructed sequences.

    Returns:
        str: The concatenated string.
    """
    res = "".join(seq +"%" for seq in spss)[:-1]
    return res + "$"