import click
import numpy

from sidh.common import attrdict
from sidh.fp import printl
from sidh.constants import strategy_data
from sympy import floor, sqrt
from math import pi


@click.command()
@click.pass_context
def csidh_parameters(ctx):
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    A = algo.curve.A
    global_L = L = algo.params.L
    n = algo.params.n
    m = algo.params.m
    delta = algo.params.delta
    full_torsion_points = algo.curve.full_torsion_points
    set_parameters_velu = algo.formula.set_parameters_velu
    set_zero_ops = algo.fp.set_zero_ops
    get_ops = algo.fp.get_ops
    KPs = algo.formula.KPs
    xISOG = algo.formula.xISOG
    xEVAL = algo.formula.xEVAL
    sJ_list = algo.formula.sJ_list
    xMUL = algo.curve.xMUL

    if True:

        # T_p belongs to E[pi - 1]
        # T_m belongs to E[pi + 1]
        T_p, T_m = full_torsion_points(A)

    else:

        # T_m belongs to E[pi - 1]
        # T_p belongs to E[pi + 1]
        T_m, T_p = full_torsion_points(A)

    assert len(L) == n
    parameters = dict()
    for idx in range(0, n, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(0, n, 1):
            if i != idx:
                Tp = xMUL(Tp, A, i)

        # -------------------------------------------------------------
        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
        # These paramters are required in KPs, xISOG, and xEVAL
        if global_L[idx] == 3:
            b = 0
            c = 0
        else:
            b = int(floor(sqrt(global_L[idx] - 1) / 2.0))
            c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))

        b += 1
        parameters[str(idx)] = []
        for j in range(0, int(floor(sqrt(pi * (b - 1)) / 1.0)), 1):

            b -= 1
            c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))
            set_parameters_velu(b, c, idx)

            total_cost = [0, 0, 0]
            # -------------------------------------------------------------
            # KPs procedure
            set_zero_ops()
            KPs(Tp, A, idx)
            t = get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # -------------------------------------------------------------
            # xISOG
            set_zero_ops()
            B = xISOG(A, idx)
            t = get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # -------------------------------------------------------------
            # xEVAL bench
            set_zero_ops()
            Tm = xEVAL(T_m, A)

            t = get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # assert(validate(B))
            parameters[str(idx)].append((b, c, t, total_cost[0] + total_cost[1]))

        if global_L[idx] == 3:
            parameters[str(idx)] = (0, 0, None, None)
        else:
            parameters[str(idx)] = min(
                parameters[str(idx)], key=lambda tup: tup[3]
            )

        print(parameters[str(idx)][0], parameters[str(idx)][1])
    return attrdict(name='parameters', **locals())

