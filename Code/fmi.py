#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions for dealing with FM-index approach for indexing a genome (coded on A,C,G,T, %)
"""
import pickle    #to serialize and deserialize objects for saving FM indexes to files.
import tools_karkkainen_sanders as tks #to sort the suffixes of the sequence


def load_fm_index(filename):
    """
    Load the FM-index from a file and return it
    """
    try:
        with open(filename, "rb") as file:
            return pickle.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File {filename} not found.")


class FmIndex:
    """
    Class for storing a FM-index of a sequence
    """

    def __init__(self, sequence):
        """
        Constructor
        """
        self.sequence = sequence
        self.sa = tks.simple_kark_sort(sequence)
        self.bwt = self.set_bwt()
        self.n, self.rank = self.set_n_and_rank()
        self.lf("T", 2)

    def save(self, filename):
        """
        Save the FM-index in a file
        """
        with open(filename, "wb") as file:
            pickle.dump(self, file)

    def set_bwt(self):
        """
        Given a sequence s and its suffix array sa, provides the associated borrows wheeler transform
        """
        bwt_list = []
        for i in self.sa:
            if i == 0:
                bwt_list.append(self.sequence[-1])
            else:
                bwt_list.append(self.sequence[i - 1])
        bwt= "".join(bwt_list)
        # print(f"BWT: {bwt}")
        # print(f"SA: {self.sa}")
        return bwt


    def set_n_and_rank(self):
        """
        N is the number of occurrences of each letter in the sequence (and so in its bwt)
        For each character in bwt, rank is the cumulated number of this caracter for each position in the bwt.
        This function computes and returns N and rank
        """
        alphabet = ["$","A", "C", "G", "T","%"]
        count = {char: 1 for char in alphabet}
        rank = []
        for char in self.bwt:
            rank.append(count[char])
            count[char] += 1
        
        n = {"A":0,"T":0,"C":0,"G":0,"$":0,"%":0}
        for l in self.bwt:
            n[l] +=1
        n = {"$":0,
             "%":1,
             "A":n["%"]+1,
             "C":n["%"]+n["A"]+1,
             "G":n["%"]+n["A"]+n["C"]+1,
             "T":n["%"]+n["A"]+n["C"]+n["G"]+1} 
        # print(f"n: {n}")
        # print(f"rank: {rank}")
        return n, rank

    def lf(self, alpha, k) -> int:
        """
        Returns the index in the suffix array corresponding to the k th suffix starting with letter alpha.
        """  
        start_pos = self.n[alpha]
        # print(f"start_pos: {start_pos}")
        return  start_pos + k - 1
        # print(f"lf: {lf}")        


    #Find the first position in the BWT that is greater than or equal to l such that BWT[position] = alpha   
    def find_next(self, alpha, l_start,l_stop) -> int:
        """
        Find the first line >= l such that BWT[line]==alpha
        """
        i = l_start 
        while i <= l_stop:
            # print(f"i: {i},r:{r}")
            if self.bwt[i] == alpha:
                return i
            i += 1
        return -1

    

    def find_prev(self, alpha, l_start,l_stop) -> int:
        """
        Find the last line <= l such that BWT[line]==alpha
        """
        i = l_stop
      
        while i >= l_start:
            #print(f"i: {i}")
            if self.bwt[i] == alpha:
                return i
            i -= 1
        return -1


    def contains(self, q) -> bool:
        """
        Check if the query q is indexed in the FM-index
        """
        l_start = 0
        l_stop = len(self.bwt) - 1
        i_pattern = len(q) - 1
        
        while i_pattern >= 0:
            # print(f"i_pattern: {i_pattern},q_pattern:{q[i_pattern]},l_start:{l_start},l_stop:{l_stop}")
            current_chr = q[i_pattern]
            next_start = self.find_next(current_chr, l_start,l_stop)
            if next_start == -1:
                return False
            next_stop = self.find_prev(current_chr, l_start,l_stop)   
            
            # print(f"2next_start: {next_start}, 2next_stop: {next_stop}")
            i_pattern -= 1
            
            l_start = self.lf(current_chr, self.rank[next_start])
            l_stop = self.lf(current_chr, self.rank[next_stop])
            # print(f"3l_start: {l_start}, 3l_stop: {l_stop}")
        return True


def main():
    """
    Main function, simple tests
    """

    my_fm_index = FmIndex("GCAC%TAGCGCAGCTAC%T$")
    # get posititions of the query "ACT"
    if my_fm_index.contains("ACT"):
        print("ACT is in the sequence")
    else:
        print("ACT is not in the sequence")
    my_fm_index.save("my_fm_index")


if __name__ == "__main__":
    main()
