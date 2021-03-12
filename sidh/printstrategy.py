"""
-----------------------------------------------------------------------------
Strategies written as a graph (a triangle a la SIDH)
-----------------------------------------------------------------------------
"""
import click
import matplotlib.pyplot as plt
import networkx as nx

from sidh.common import attrdict, strategy_evaluation
from sidh.constants import strategy_data

@click.command()
@click.pass_context
def print_strategy(ctx):
    """
    draw graphs
    """
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    if setting.algorithm == 'csidh':
        n = algo.params.n
        m = algo.params.m
        L = algo.params.L
        rounds = algo.gae.rounds
        geometric_serie = algo.gae.geometric_serie

        BATCHES = lambda n, sigma: [
            len([j for j in range(i, n, sigma)]) for i in range(sigma)
        ]


        if len(set(m)) > 1:
            # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
            LABEL_m = 'different_bounds'
        else:
            # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
            LABEL_m = 'with_same_bounds'

        if setting.verbose:
            verb = '-suitable'
        else:
            verb = '-classical'

        # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
        m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
        r_out, L_out, R_out = rounds(m_prime[::-1], n)
        for j in range(0, len(r_out), 1):

            R_out[j] = list([L[::-1][k] for k in R_out[j]])
            L_out[j] = list([L[::-1][k] for k in L_out[j]])

        f = open(
            strategy_data
            + setting.algorithm
            + '-'
            + setting.prime
            + '-'
            + setting.style
            + '-'
            + setting.formula
            + '-'
            + LABEL_m
            + verb
        )
        print("// Strategies to be read from a file")
        S_out = []
        for i in range(0, len(r_out), 1):

            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            S_out.append(tmp)

        f.close()

        filename = ( setting.algorithm + '-' + setting.prime + '-' + setting.style
                     + '-' + setting.formula + '-' + LABEL_m + verb)

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
                edge_labels=True,
            )

            # Saving the graph as a .PNG figure
            plt.savefig(filename + '-id' + str(idx) + '.png')
            print("saving graph: " + filename + '-id' + str(idx) + '.png')
            # plt.show()
            plt.close()

        print("// The strategies have been plotted and stored in the current directory")
    else:
        print("bsidh not yet implemented")
        click.Exit(1)

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

    def simba(n, sigma):

        vertexes = dict()  # list of the position of each node
        edges = []  # edges of the isogeny triangle

        edge_colors = (
            []
        )  # color of each edge: blue for scalar multiplications, and orange for isogeny evalutions
        vertex_colors = (
            []
        )  # color of each node: red for the leaves, otherwise color is set to white

        batches = BATCHES(n, sigma)
        t = -1
        acc = 0
        for i in range(sigma):

            t += 1
            vertexes[t] = (acc, 0)  # Initial vertex of the current batch
            vertex_colors.append('black')

            t += 1
            vertexes[t] = (acc, -n + batches[i])  # Root of the current batch
            vertex_colors.append('black')

            edges.append((t - 1, t))  # New edge to be added
            edge_colors.append(
                'tab:blue'
            )  # Color of this type of edge is always blue

            vertex = vertexes[t]
            for j in range(batches[i] - 1):

                t += 1
                vertexes[t] = (
                    vertexes[t - 1][0],
                    vertexes[t - 1][1] - batches[i] + j + 1,
                )  # Next leaf
                vertex_colors.append('black')

                edges.append((t - 1, t))  # New edge to be added
                edge_colors.append(
                    'tab:blue'
                )  # Color of this type of edge is always blue

                t += 1
                vertexes[t] = (vertexes[t - 2][0] + 1, vertexes[t - 2][1])
                vertex_colors.append('black')

                edges.append((t - 2, t))  # New edge to be added
                edge_colors.append(
                    'tab:red'
                )  # Color of this type of edge is always blue

            acc += batches[i]

        return vertexes, vertex_colors, edges, edge_colors
    return attrdict(name='print-strategy', **locals())
