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
import statistics
import argparse
import os
import sys
import scipy
import random
from random import randint
import networkx as nx
# import matplotlib
# import matplotlib.pyplot as plt
from operator import itemgetter
random.seed(9001)



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
    """
    Parameters:
    -----------
    seq : string
    -> sequence a découpé
    Return:
    -------
    k : integer
    -> taille du kmer
    """
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]


def build_kmer_dict(fichier : str, k : int):
    """
    Parameters:
    -----------
    fichier : str
    -> chemin du fichier ,
    K : int
    -> taille du kmers
    Return:
    -------
    out : dictionnary
    -> dictionaire de Kmers
    """
    out = {}
    for i in read_fastq(fichier):
        for j in cut_kmer(i,k = k ):
            if not j in out :
                out.setdefault(j, 1)
            else :
                out[j] += 1
    return out


def build_graph(dic):
    """
    Parameters:
    -----------
    dic : dictionaire de kmers
    -> dictionaire des kmers
    Return:
    -------
    graph : nx.DiGraph
    -> graphe orienté construit a partir de dic
    """
    graph = nx.DiGraph()
    for i in dic:
        start = i[1:]
        end = i[:-1]
        graph.add_node(start)
        graph.add_node(end)
        graph.add_edge(end,start,weight = dic[i])
    return graph


def get_starting_nodes(graph : nx.DiGraph):
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph des kmers
    Return:
    -------
    start : list
    -> list des noeud d'entré
    """
    start = []
    for i in graph.nodes:
        pred = graph.predecessors(i)
        if len(list(pred)) == 0 :
            start.append(i)
    return start


def get_sink_nodes(graph : nx.DiGraph):
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph des kmers
    Return:
    -------
    out : list
    -> list des noeud de sortie
    """
    out = []
    for i in graph.nodes:
        succ = list(graph.successors(i))
        if len(succ) == 0 :
            out.append(i)
    return out


def get_contigs(graph, in_nodes , out_nodes):
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graphe de kmers
    in_nodes : list
    -> liste des noeuds d'entre
    out_nodes : list
    -> liste des noeud de sortie
    Return:
    -------
    contigs : list
    -> list de tuple (contigs, taille du contigs)

    """
    contigs = []
    for in_ in in_nodes:
        for out_ in out_nodes:
            for path in nx.all_simple_paths(graph, in_, out_):
                seq = ""
                seq += path[0]
                for i in path[1:]:
                    seq += i[-1]
                contigs.append((seq, len(seq) ))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs, out_file):
    """
    Parameters:
    -----------
    contigs : list
    -> list de tuple(contigs , taille des contigs)
    out_file : string
    -> chemin jusqu'au fichier

    Return:
    None
    -------

    """
    with open(out_file, "w") as filout:
        for i,val in enumerate(contigs):
            filout.write(">contig_{} len={}\n".format(i, val[1]))
            filout.write("{}\n".format(fill(val[0])))


def std(list_in):
    """
    Parameters:
    -----------
    list_in : list
    -> list de valeur
    Return:
    -------
    std : double
    -> standant variaiton de list_in

    """
    return statistics.stdev(list_in)


def path_average_weight(graph , path):
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graphe des kmers
    path : list
    -> list des noeud composant le chemin a élavuer

    Return:
    -------
    mean : double
    -> moyenne des poids du chemin du graph

    """
    print("Path:", path)
    vals = []
    for i,val  in enumerate(path[:-1]):
        edge_d= graph.get_edge_data(val, path[i+1])
        vals.append(edge_d["weight"])
    return statistics.mean(vals)


def remove_paths(graph, paths,
                 delete_entry_node: bool,
                 delete_sink_node: bool):
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph de kmers
    paths : list
    -> list des chemins a retiré
    delete_entry_node : bool
    -> specifiier si le premier noeud du chemnin doit etre retiré
    delete_sink_node : bool
    -> spécifier si le dernier noued du chemein doit etre détruit

    Return:
    -------
    grpah : nx.graph
    -> graphe nettoyer des noeuds des chemins

    """
    start = 0 if delete_entry_node else 1
    end = 0 if delete_sink_node else 1
    for i, val  in enumerate(paths):
        for j in val[start:len(val)-end]:
            graph.remove_node(j)
    return graph



def select_best_path(graph, paths, length, weight,
                     delete_entry_node = False,
                     delete_sink_node = False):
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph de kmers
    paths : list
    -> list des chemins a retiré
    length : list
    -> list des taille associer au graphe de paths
    weight : list
    -> list des poid associer au chemin de paths
    delete_entry_node : bool
    -> specifiier si le premier noeud du chemnin doit etre retiré
    delete_sink_node : bool
    -> spécifier si le dernier noued du chemein doit etre détruit

    Return:
    -------
    grpah : nx.graph
    -> graphe nettoyer des noeuds des chemin les moins interresant
    """
    if weight[0] > weight[1]:
        graph = remove_paths(graph,[paths[1]],delete_entry_node,delete_sink_node)
    elif weight[0]< weight[1]:
        graph = remove_paths(graph,[paths[0]], delete_entry_node,delete_sink_node)
    elif weight[0] == weight[1]:
        if length[0] > length[1]:
            graph = remove_paths(graph,[paths[1]],delete_entry_node,delete_sink_node)
        elif length[0] < length[1]:
            graph = remove_paths(graph,[paths[0]], delete_entry_node,delete_sink_node)
        elif length[0] == length[1]:
            rand = random.randint(0,1)
            graph = remove_paths(graph,[paths[rand]], delete_entry_node,delete_sink_node)
    return graph



def solve_bubble(graph,start,end) :
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph de kmers
    start : node
    -> noeud de début de la bulle
    end : node
    -> node de fin de la bulle
    Return:
    -------
    grpah : nx.graph
    -> graph nettoyer de la bulle
    """
    while len(list(nx.all_simple_paths(graph,start,end))) > 1 :
        paths = list(nx.all_simple_paths(graph,start,end))
        mean = [path_average_weight(graph, paths[0]),
                path_average_weight(graph, paths[1])]
        lengths = [len(paths[0]),len(paths[1])]
        graph = select_best_path(graph, paths[0:2],lengths,mean)
    return graph


def get_buble(path_1, path_2):
    """
    Parameters:
    -----------
    path_1 : list
    -> chemin d'une bulle
    path_2 : lsit
    -> chemin 2 de la bulle

    Return:
    -------
    p1 : node
    -> debut d'une sous bulle
    p2 : node
    -> node de fin d'une sous bulle
    """
    pos_1=None
    for i,node in enumerate(path_1):
        if node not in path_2:
            pos_1 = i -1
            break
    pos_2 = None
    for i , node in enumerate(path_1[pos_1+2:]):
        if node in path_2:
            print("ok",node)
            pos_2 = i + pos_1 + 2
            break
    if pos_1 == None or pos_2 == None:
        return None, None
    return path_1[pos_1],path_1[pos_2]



def find_sub_buble(graph, node_in , node_out,path_1, path_2):
    print("activate")
    cond = True
    print("out", node_out)

    for i in graph.predecessors(node_out):
        if i in path_1 or i in path_2 : 
            print(i , list(nx.all_simple_paths(graph, node_in, i)))
            if len(list(nx.all_simple_paths(graph, node_in, i))) != 1 :
                print("stop")
                cond = False
    if cond :
        print("ok")
        return(node_in, node_out)
    else:
        print("else")
        node_1 = list(graph.successors(node_in))[0]
        node_2 = list(graph.predecessors(node_out))[0]
        paths = list(nx.all_simple_paths(graph,node_in,node_2))
        if len(paths) == 2:
            return (node_in, node_2)
        elif len(paths) > 2:
            return find_sub_buble(graph,node_in,node_2,paths[0],paths[1])
        paths = list(nx.all_simple_paths(graph,node_1,node_out))
        if len(paths) == 2:
            return (node_1, node_out)
        elif len(paths) > 2:
            return find_sub_buble(graph,node_1, node_out, paths[0], paths[1])


def simplify_bubbles(graph) :
    """
    Parameters:
    -----------
    graph : nx.graph
    -> graph des kmers

    Return:
    -------
    graph : nx.graph
    -> graphe nétoyer de bulle
    """
    print()
    in_ = get_starting_nodes(graph)
    out_ = get_sink_nodes(graph)
    for i in in_:
        for j in out_:
            while len(list(nx.all_simple_paths(graph,i,j))) > 1:
                paths = list(nx.all_simple_paths(graph,i,j))
                path_1 = paths[0]
                path_2 = paths[1]
                print("#############################")
                print(path_1)
                print(path_2)
                node_1 = path_1[0]
                node_2 = path_2[-1]
                n_1, n_2 = get_buble(path_1[path_1.index(node_1):path_1.index(node_2)+1],
                                     path_2[path_1.index(node_1):path_1.index(node_2)+1])
                # print((n_1,n_2))
                # #n_1,n_2 = find_sub_buble(graph, n_1, n_2)
                # n_1, n_2 = find_sub_buble(graph, n_1, n_2,path_1,path_2)
                graph = solve_bubble(graph, n_1,n_2)
    return graph

def solve_entry_tips(graph ,entre) :
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph des kmers
    entre : list
    -> list des noeud d'entre
    Return:
    -------
    graph : nx.DiGraph
    -> graphe nettoyé des neoud d'entre
    """

    while len(get_starting_nodes(graph)) > 1:
        starts = get_starting_nodes(graph)
        start = starts[0]
        node = start
        while True:
            if len(list(graph.predecessors(node))) > 1:
                break
            node = list(graph.successors(node))[0]
        n_1 = list(graph.predecessors(node))[0]
        n_2 = list(graph.predecessors(node))[1]
        explor_1 = True
        explor_2 = True
        to_explore= n_1
        join = node
        while True :
            if not explor_1 and not  explor_2:
                break
            if explor_1:
                pred = list(graph.predecessors(n_1))
                if len(pred) ==1 :
                    n_1 = pred[0]
                elif len(pred) == 0 :
                    explor_1 = False
                else:
                    node = n_1
                    n_1 = pred[0]
                    n_2 = pred[1]
            elif explor_2 :
                pred = list(graph.predecessors(n_2))
                if len(pred) == 1:
                    n_2 = pred[0]
                elif len(pred) == 0 :
                    explor_2 = False
                else :
                    node = n_2
                    n_1 = pred[0]
                    explor_1 = True
                    n_2 = pred[1]
        path_1 = list(nx.all_simple_paths(graph,n_1,join))[0]
        path_2 = list(nx.all_simple_paths(graph,n_2,join))[0]
        w_1 = path_average_weight(graph,path_1)
        w_2 = path_average_weight(graph,path_2)
        graph = select_best_path(graph, [path_1, path_2],
                                 [len(path_1),len(path_2)],
                                 [w_1, w_2],True,False)

    return graph


def solve_out_tips(graph, out) :
    """
    Parameters:
    -----------
    graph : nx.DiGraph
    -> graph des kmers
    entre : list
    -> list des noeud de sortie
    Return:
    -------
    graph : nx.DiGraph
    -> graphe nettoyé des neoud de sortie
    """

    while len(get_sink_nodes(graph)) > 1:
        ends = get_sink_nodes(graph)
        end = ends[0]
        node = end
        while True:
            if len(list(graph.successors(node))) > 1:
                break
            node = list(graph.predecessors(node))[0]
        n_1 = list(graph.successors(node))[0]
        n_2 = list(graph.successors(node))[1]
        explor_1 = True
        explor_2 = True
        to_explore= n_1
        join = node
        while True :
            if not explor_1 and not  explor_2:
                break
            if explor_1:
                pred = list(graph.successors(n_1))
                if len(pred) ==1 :
                    n_1 = pred[0]
                elif len(pred) == 0 :
                    explor_1 = False
                else:
                    node = n_1
                    n_1 = pred[0]
                    n_2 = pred[1]
            elif explor_2 :
                pred = list(graph.successors(n_2))
                if len(pred) == 1:
                    n_2 = pred[0]
                elif len(pred) == 0 :
                    explor_2 = False
                else :
                    node = n_2
                    n_1 = pred[0]
                    explor_1 = True
                    n_2 = pred[1]
        path_1 = list(nx.all_simple_paths(graph,join,n_1,))[0]
        path_2 = list(nx.all_simple_paths(graph,join,n_2))[0]
        w_1 = path_average_weight(graph,path_1)
        w_2 = path_average_weight(graph,path_2)
        graph = select_best_path(graph, [path_1, path_2],
                                 [len(path_1),len(path_2)],
                                 [w_1, w_2],False,True)
    
    return graph


def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    k_dict = build_kmer_dict(args.fastq_file,
                             k =args.kmer_size)
    print("k_dict")
    graph = build_graph(k_dict)
    print("graph")
    starts = get_starting_nodes(graph)
    ends = get_sink_nodes(graph)
    print("end","start")
    graph = simplify_bubbles(graph)
    for i, start in enumerate(starts[:-1]):
        graph = solve_entry_tips(graph,[start,starts[i+1]])
    for i , end in enumerate(ends[:-1]):
        graph = solve_out_tips(graph, [end,ends[i+1]])

    starts = get_starting_nodes(graph)
    ends = get_sink_nodes(graph)
    contigs = get_contigs(graph, starts, ends)
    save_contigs(contigs, args.output_file )



if __name__ == '__main__':
    main()

    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from([(3, 2, 10), (2, 4, 15), (4, 5, 15),
                                     (2, 10,10), (10, 5,10), (2, 8, 3),
                                     (8, 9, 3), (9, 5, 3), (5, 6, 10),
                                     (5, 7, 10)])
    graph_1 = simplify_bubbles(graph_1)

    print( (2,8) not in graph_1.edges())
    print( (8,9) not in graph_1.edges())
    print( (9,5) not in graph_1.edges())
    print( (2,10) not in graph_1.edges())
    print( (10, 5) not in graph_1.edges())
