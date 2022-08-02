from functools import wraps
from math import floor

# Dictionary which provides attribute access to its keys.
class attrdict(dict):
    __getattr__ = dict.__getitem__

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
        if self.__class__ is not other.__class__:
            other = self.__class__(other)
        else:
            if self.__class__.p != other.__class__.p:
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

def geometric_serie(m, l):
    """
    geometric_serie()
    inputs: and integer m, and a prime number l
    output: the nearest integer to
                  l
            m x -----
                l - 1
    """
    l_float = float(l)
    m_float = float(m)
    return floor((m_float * l_float) / (l_float - 1.0) + 0.5)

def filtered(List, sublist):
    """
    filtered()
    inputs : a list L and a sublist SL of L
    output : L \\ SL
    """
    return [e for e in List if e not in sublist]

def rounds(e, n):
    """
    rounds()
    inputs : an integer vector (maximum number of isogeny constructions to be performed),
             and the length of the vector
    output : the subset of (indexes of) the small odd primes that determines the optimal
             strategy to be used, the number of times that each strategy will be used, and
             the complement of each subset (with respect to the set of all the small odd primes)
    """
    tmp_N = range(n)
    tmp_e = list(e)
    rounds_out = []
    sublists_L = []
    sublists_C = []
    while [e_i for e_i in tmp_e if e_i > 0] != []:
        e_min = min([e_i for e_i in tmp_e if e_i > 0])
        rounds_out.append(e_min)
        sublists_L.append([i for i in tmp_N if tmp_e[i] >= e_min])
        sublists_C.append(filtered(tmp_N, sublists_L[len(sublists_L) - 1]))
        tmp_e = [(tmp_e[i] - e_min) for i in tmp_N]

    return rounds_out, sublists_L, sublists_C
