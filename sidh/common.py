class attrdict(dict):
    """
    Dictionary which provides attribute access to its keys.
    """

    def __getattr__(self, key):
        if key in self:
            return self[key]
        else:
            raise AttributeError(
                "%r object has no attribute %r" % (type(self).__name__, key)
            )

def strategy_evaluation(strategy, n):

    vertexes = dict()  # list of the position of each node
    edges = []  # edges of the isogeny triangle

    edge_colors = (
        []
    )  # color of each edge: blue for scalar multiplications, and orange for isogeny evalutions
    vertex_colors = (
        []
    )  # color of each node: red for the leaves, otherwise color is set to white

    vertexes[0] = (0.0, 0.0)  # Root of the strategy
    ramifications = []  # nodes having two different types of edges
    moves = [0]  #
    k = 0  # current element of the strategy
    t = 0  # current vertex

    # Strategy evaluation starts next
    ramifications.append([0, vertexes[0]])  # The root is always a ramification
    for i in range(len(strategy)):

        # Getting the new possible ramifications
        while sum(moves) < (n - 1 - i):

            vertex = ramifications[-1][1]

            moves.append(strategy[k])  # Increasing moves
            k += 1  # Next element of the strategy

            t += 1  # Next vertex
            edges.append((t - 1, t))  # New edge to be added
            edge_colors.append(
                'tab:blue'
            )  # Color of this type of edge is always blue

            vertexes[t] = (
                i,
                vertex[1] - strategy[k - 1],
            )  # New vertex to be added
            ramifications.append([t, vertexes[t]])  # New ramification
            vertex_colors.append('black')

        # Removing the ramifications (it is not required more!)
        ramifications.pop()
        moves.pop()

        # Updating the ramifications
        for j in range(len(ramifications)):

            t += 1
            vertexes[t] = (
                ramifications[j][1][0] + 1.0,
                ramifications[j][1][1],
            )
            edges.append((ramifications[j][0], t))

            ramifications[j] = [t, vertexes[t]]
            edge_colors.append('tab:red')

    return vertexes, vertex_colors, edges, edge_colors
