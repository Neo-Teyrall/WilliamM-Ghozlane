#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from operator import itemgetter
import random
import scipy
random.seed(9001)
from random import randint
import statistics

__author__ = "William margerit"
__copyright__ = "Universite Paris"
__credits__ = ["William Margerit"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "William Margerit"
__email__ = "mailliw.marg@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


#==============================================================
# Main program
#==============================================================

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

def build_graph(dic):
    G = nx.DiGraph()
    for i in dic:
        s = i[1:]
        p = i[:-1]
        G.add_node(s)
        G.add_node(p)
        G.add_edge(p,s,weight = dic[i])
    return G

def get_starting_nodes(G : nx.DiGraph):
    start = []
    for i in G.nodes:
        a = G.predecessors(i)
        if len(list(a)) == 0 :
            start.append(i)
    return start

def get_sink_nodes(G : nx.DiGraph):
    out = []
    for i in G.nodes:
        a = list(G.successors(i))
        if len(a) == 0 :
            out.append(i)
    return out

def get_contigs(G, in_nodes , out_nodes):
    contigs = []
    for in_ in in_nodes:
        for out_ in out_nodes:
            for path in nx.all_simple_paths(G, in_, out_):
                contigs.append((path, len(path)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs, out_file):
    with open(out_file, "w") as filout:
        for i,val in enumerate(contigs):
            filout.write(">contigs_{} len={}\n".format(i, val[1]))
            seq = val[0][0]
            for j in val[0][1:]:
                seq += j[-1]
            #seq = fill(seq)
            filout.write("{}\n".format(seq))
    pass

def main():
    """
    Main program function
    """
    # Get arguments
    #args = get_arguments()
    
    k_dict = build_kmer_dict("../data/eva71_two_reads.fq", k = 6)
    G = build_graph(k_dict)
    start_point = get_starting_nodes(G)
    end_nodes = get_sink_nodes(G)
    contigs = get_contigs(G, start_point, end_nodes)
    save_contigs(contigs, "out.txt" )

def save_contigs() :
    pass 
def std() : 
    pass
def path_average_weight() : 
    pass
def remove_paths() : 
    pass
def delete_entry_node() :
    pass
def select_best_path() : 
    pass
def solve_bubble() :
    pass
def simplify_bubbles() : 
    pass
def solve_entry_tips() : 
    pass
def solve_out_tips() : 
    pass
    

if __name__ == '__main__':
    main()
