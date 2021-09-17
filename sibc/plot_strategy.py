"""
-----------------------------------------------------------------------------
Strategies written as a graph (a triangle a la SIDH)
-----------------------------------------------------------------------------
"""
import sys
import click
import matplotlib.pyplot as plt
import networkx as nx

from sibc.common import attrdict, strategy_evaluation, geometric_serie, rounds
from pkg_resources import resource_filename

@click.command()
@click.pass_context
def plot_strategy(ctx):
    """
    draw strategy graphs as a subgraph Discrete Right Triangle (DRT)
    """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    if setting.algorithm == 'csidh':
        L = algo.params.L
        m = algo.params.m
        n = algo.params.n
        
        temporal_m = list(set(m))
        if len(temporal_m) > 1:
            # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
            bounds = '-diffbounds'
        else:
            # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
            bounds = '-samebounds'

        # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
        m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
        r_out, L_out, R_out = rounds(m_prime[::-1], n)
        for j in range(0, len(r_out), 1):

            R_out[j] = list([L[::-1][k] for k in R_out[j]])
            L_out[j] = list([L[::-1][k] for k in L_out[j]])

        file_path = (
            "data/strategies/"
            + algo.curve.model
            + '/csidh/'
            + 'csidh'
            + '-'
            + setting.prime
            + '-'
            + setting.style
            + '-e'
            + setting.exponent
            + bounds
            + '-'
            + setting.formula
            + '-'
            + algo.formula.multievaluation_name
            + algo.formula.tuned_name
        )
        file_path = resource_filename('sibc', file_path)
        f = open(file_path)
        S_out = []
        for i in range(0, len(r_out), 1):

            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            S_out.append(tmp)

        f.close()

        file_path = (
            'csidh'
            + '-'
            + setting.prime
            + '-'
            + setting.style
            + '-e'
            + setting.exponent
            + bounds
            + '-'
            + setting.formula
            + '-'
            + algo.formula.multievaluation_name
            + algo.formula.tuned_name
            + '-'
            + algo.curve.model
        )

    elif setting.algorithm == 'bsidh':
        file_path = (
            "data/strategies/"
            + algo.curve.model
            + '/bsidh/'
            + 'bsidh'
            + '-'
            + setting.prime
            + '-'
            + setting.formula
            + '-'
            + algo.formula.multievaluation_name
            + algo.formula.tuned_name
        )
        file_path = resource_filename('sibc', file_path)
        f = open(file_path)
        S_out = []
        # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ...,
        # l_{n-1}]
        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        S_out.append(list(tmp))
        # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ...,
        # l_{n-1}]
        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        S_out.append(list(tmp))
        f.close()

        file_path = (
            'bsidh'
            + '-'
            + setting.prime
            + '-'
            + setting.style
            + '-'
            + setting.formula
            + '-'
            + algo.formula.multievaluation_name
            + algo.formula.tuned_name
            + '-'
            + algo.curve.model
        )

    elif setting.algorithm == 'sidh':
        file_path = (
                "data/strategies/"
                + algo.curve.model
                + '/sidh/'
                + 'sidh'
                + '-'
                + setting.prime
        )
        file_path = resource_filename('sibc', file_path)
        f = open(file_path)
        S_out = []
        # Corresponding to 2-isogenies
        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        S_out.append(list(tmp))
        # Corresponding to 3-isogenies
        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        S_out.append(list(tmp))
        f.close()

        file_path = (
                'sidh'
                + '-'
                + setting.prime
                + '-'
                + algo.curve.model
        )

    else:
        print("only csidh, bsidh, and sidh are implemented")
        sys.exit(1)

    # ----
    for idx in range(0, len(S_out), 1):
        S = S_out[idx]
        n = len(S) + 1

        # Strategy written as a graph
        vertexes, vertex_colors, edges, edge_colors = strategy_evaluation(S, n)

        # Simba method written as a graph
        # vertexes, vertex_colors, edges, edge_colors = simba(n, 3)

        # All the Discrete Right Triangle
        # vertexes, vertex_colors, edges, edge_colors = DRT(n)
        G = nx.Graph()

        # Adding nodes in specific positions
        G.add_nodes_from(list(range(len(vertexes))))

        nx.set_node_attributes(G, vertexes, 'pos')
        # Adding edges with specific colors
        for i in range(len(edges)):
            G.add_edge(edges[i][0], edges[i][1], color=edge_colors[i])

        # Setting variables for a pretty plot of the graph
        edges = G.edges()
        edge_colors = [G[u][v]['color'] for u, v in edges]
        weights = [6 for u, v in edges]
        vertex_sizes = [24] * len(vertexes)

        # Finally, the graph will be plotted
        plt.figure(1, figsize=(17, 17))
        nx.draw(
            G,
            vertexes,
            node_color=['black'] * len(vertexes),
            node_size=vertex_sizes,
            edge_color=edge_colors,
            width=weights,
            #edge_labels=True,
        )

        # Saving the graph as a .PNG figure
        file_name = file_path + '-id' + str(idx) + '.png'
        plt.savefig(file_name)
        print("saving graph: " + file_name)
        # plt.show()
        plt.close()

    print("// The strategies have been plotted and stored in the current directory")

    def DRT(n):

        vertexes = dict()  # list of the position of each node
        vertex_colors = (
            []
        )  # color of each node: red for the leaves, otherwise color is set to white
        acc = 0

        # Different shape of the isogeny graph
        for i in range(n):
            for j in range(n - 1 - i):
                vertex_colors.append('black')
                vertexes[acc] = (i, -j)
                acc += 1

        return vertexes, vertex_colors, [], []

    return attrdict(name='plot-strategy', **locals())
