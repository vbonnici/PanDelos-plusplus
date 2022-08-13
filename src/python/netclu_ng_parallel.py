#!/usr/bin/python3

import sys
import networkx as nx
from multiprocessing import Pool
import itertools

def girvan_newman(G, most_valuable_edge=None):

    if G.number_of_edges() == 0:
        yield tuple(nx.connected_components(G))
        return
    # If no function is provided for computing the most valuable edge,
    # use the edge betweenness centrality.
    if most_valuable_edge is None:

        def chunks(l, n):
            """Divide a list of nodes `l` in `n` chunks"""
            l_c = iter(l)
            while 1:
                x = tuple(itertools.islice(l_c, n))
                if not x:
                    return
                yield x

        def most_valuable_edge(G):
            """Returns the edge with the highest betweenness centrality
            in the graph `G`.

            We have guaranteed that the graph is non-empty, so this dictionary will never be empty.

            """


            print("Parallel betweenness centrality")
            p = Pool(processes=None) # aggiungere numero processi dinamico
            node_divisor = len(p._pool) * 4
            node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
            num_chunks = len(node_chunks)

            bt_sc = p.starmap(
                nx.edge_betweenness_centrality_subset,
                zip(
                    [G] * num_chunks,
                    node_chunks,
                    [list(G)] * num_chunks,
                    [True] * num_chunks,
                    ['weight'] * num_chunks,
                    ),
            )

            betweenness = bt_sc[0]
            for bt in bt_sc[1:]:
                for n in bt:
                    betweenness[n] += bt[n]

            return max(betweenness, key=betweenness.get)

    # The copy of G here must include the edge weight data.
    g = G.copy().to_undirected()
    # Self-loops must be removed because their removal has no effect on
    # the connected components of the graph.
    g.remove_edges_from(nx.selfloop_edges(g))
    while g.number_of_edges() > 0:
        yield _without_most_central_edges(g, most_valuable_edge)


def _without_most_central_edges(G, most_valuable_edge):
    """Returns the connected components of the graph that results from
    repeatedly removing the most "valuable" edge in the graph.

    `G` must be a non-empty graph. This function modifies the graph `G`
    in-place; that is, it removes edges on the graph `G`.

    `most_valuable_edge` is a function that takes the graph `G` as input
    (or a subgraph with one or more edges of `G` removed) and returns an
    edge. That edge will be removed and this process will be repeated
    until the number of connected components in the graph increases.

    """
    original_num_components = nx.number_connected_components(G)
    num_new_components = original_num_components
    while num_new_components <= original_num_components:
        edge = most_valuable_edge(G)
        G.remove_edge(*edge)
        new_components = tuple(nx.connected_components(G))
        num_new_components = len(new_components)
    return new_components


def get_max_collision(coco, pnet):
    print("get_max_collision")
    collisions = dict()
    for s in coco:
        if seq_genome[s] not in collisions:
            collisions[seq_genome[s]] = list()
        collisions[seq_genome[s]].append(s)
    max_k = 0
    for k, v in collisions.items():
        for s1 in v:
            s_k = 0;
            for s2 in v:
                if (s1 != s2) and not (pnet.has_edge(s1, s2) or pnet.has_edge(s2, s1)):
                    s_k += 1
            if s_k > max_k:
                max_k = s_k
    # if max_k > 0:
    #    print(max_k, len(coco))
    return max_k

def split_until_max_k(coco, pnet):
    print("split_until_max_k")
    snet = pnet.subgraph(coco)

    if (len(coco) > 100):
        gcoms = girvan_newman(snet)
    else:
        gcoms = nx.algorithms.community.centrality.girvan_newman(snet)

    # coms = tuple( c for c in next(gcoms) )
    coms = tuple(sorted(c) for c in next(gcoms))
    print("gn", coms)
    rcoms = list()
    nof_coms = 0
    for com in coms:
        if get_max_collision(com, snet) > 0:
            rcoms = rcoms + split_until_max_k(com, snet)
        else:
            rcoms.append(com)
    return rcoms


def print_family(fam):
    print(len(fam))
    print("fam", sorted(fam))
    print("F{ ", end='', sep='')
    print(" ; ".join([seq_names[f] for f in sorted(fam)]), end='', sep='')
    # for f in fam:
    #        print(seq_names[f] + " ; ", end='',sep='')
    print('}\n', end='')


def print_family_descriptions(fam):
    print("D{ ", end='', sep='')
    print(" ; ".join([seq_descr[f] for f in sorted(fam)]), end='', sep='')
    # for f in fam:
    #        print(seq_descr[f] + " ; ", end='',sep='')
    print('}\n', end='')
    print("S{ ", end='', sep='')
    print(" ; ".join([d for d in set([seq_descr[f] for f in sorted(fam)])]), end='', sep='')
    # for f in fam:
    #        print(seq_descr[f] + " ; ", end='',sep='')
    print('}\n', end='')
    # print( set( [ seq_descr[f] for f in fam ]) )
    print('-')


if __name__ == '__main__':

    iseqs = sys.argv[1]
    inet = sys.argv[2]

    # seq_names = list()
    seq_names = dict()  # id -> name
    seq_genome = dict()
    seq_descr = dict()
    genomes = dict()

    i = 0
    for line in open(iseqs, 'r'):
        if i % 2 == 0:
            cols = line.strip().split('\t')
            seq_id = len(seq_names)
            #        seq_names.append(cols[1])
            seq_names[seq_id] = cols[1]
            seq_genome[seq_id] = cols[0]
            seq_descr[seq_id] = cols[2]
            if cols[0] not in genomes:
                genomes[cols[0]] = list()
            genomes[cols[0]].append(seq_id)
        i += 1

    print("nof sequences", len(seq_names))
    print("nof genomes", len(genomes))

    ccheck = set()
    for k, v in seq_names.items():
        if v in ccheck:
            print("duplicated seq name", v)
        ccheck.add(v)

    inodes = set()
    pnet = nx.Graph()
    for line in open(inet, 'r'):
        cols = line.strip().split('\t')
        cols[0] = int(cols[0])
        cols[1] = int(cols[1])
        cols[2] = float(cols[2])
        if cols[0] not in inodes:
            inodes.add(cols[0])
            pnet.add_node(cols[0])
        if (cols[1] not in inodes) and (cols[0] != cols[1]):
            inodes.add(cols[1])
            pnet.add_node(cols[1])
        if cols[0] != cols[1]:
            pnet.add_edge(cols[0], cols[1], weight=cols[2])
            pnet.add_edge(cols[1], cols[0], weight=cols[2])

    print("nof net nodes", pnet.number_of_nodes())
    print("nof net edges", pnet.number_of_edges())

    print('-' * 40)

    comps_size_distr = dict()
    comps = nx.algorithms.components.connected_components(pnet)
    nof_comps = 0
    # print('nof connected components', len(comps))
    for comp in comps:
        comps_size_distr[len(comp)] = comps_size_distr.get(len(comp), 0) + 1
        nof_comps += 1
    for k, v in sorted(comps_size_distr.items()):
        print(k, v)
    print('nof connected components', nof_comps)

    print('-' * 40)

    remaining_singletons = set()
    for g in seq_names.keys():
        remaining_singletons.add(g)

    fnodes = set()
    coms_size_distr = dict()
    nof_coms = 0
    for coco in nx.algorithms.components.connected_components(pnet):
        print('-' * 10)
        print("coco", sorted(coco))
        max_k = get_max_collision(coco, pnet)
        if max_k > 0:
            print(max_k, len(coco))
            coms = split_until_max_k(coco, pnet)
            # print("coms", sorted(coms))
            nof_coms += len(coms)
            for com in coms:
                coms_size_distr[len(com)] = coms_size_distr.get(len(com), 0) + 1
                # print("com", sorted(com))
                print_family(com)
                for g in com:
                    remaining_singletons.remove(g)
                print_family_descriptions(com)
        else:
            nof_coms += 1
            coms_size_distr[len(coco)] = coms_size_distr.get(len(coco), 0) + 1
            # print("coco", sorted(coco))
            print_family(coco)
            for g in coco:
                remaining_singletons.remove(g)
            print_family_descriptions(coco)

    for g in remaining_singletons:
        print("F{ " + seq_names[g] + " }")

    for k, v in sorted(coms_size_distr.items()):
        print(k, v)
    print('nof communities', nof_coms)
    print('-' * 40)