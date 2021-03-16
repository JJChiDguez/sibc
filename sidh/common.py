from functools import wraps

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

def check(func):
    @wraps(func)
    def method(self, other):
        if type(self) is not type(other):
            other = self.__class__(other)
        else:
            if self.field.p != other.field.p:
                raise ValueError
        return func(self, other)

    return method

def doc(s):
    class __doc(object):
        def __init__(self,f):
            self.func = f
            self.desc = s
        def __call__(self,*args,**kwargs):
            return self.func(*args,**kwargs)
        def __repr__(self):
            return self.desc
    return __doc

'''
    chunks()
    inputs: a string, a list, and the maximum  number of elements in each chunk
    -----
    NOTE: This function divide the input list into len(L) / k chunks.
'''
chunks = (
    lambda NAME, L, n: [NAME + ' =\t{']
    + [
        '\t' + ','.join(list(map(format, L[i * n : (i + 1) * n], ['3d'] * n)))
        for i in range((len(L) + n - 1) // n)
    ]
    + ['\t};']
)
'''
    printl()
    inputs: a string, a list, and the maximum number k of elements in each chunk
    -----
    NOTE: this function prints a given list by chunks of size k.
'''


def printl(NAME, L, k):

    to_print = chunks(NAME, L, k)
    print(to_print[0])
    for i in range(1, len(to_print) - 2):
        print(to_print[i] + ",")

    print(to_print[len(to_print) - 2])
    print(to_print[len(to_print) - 1])