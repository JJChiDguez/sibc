import sys
import click
from pkg_resources import resource_filename
import numpy

from sibc.common import attrdict, printl
from math import floor, sqrt, pi


@click.command()
@click.pass_context
def csidh_precompute_parameters(ctx):
    """ Precomputation of tuned velusqrt parameters """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    L = algo.params.L
    n = algo.params.n
    m = algo.params.m
    generators = algo.curve.generators
    set_parameters_velu = algo.formula.set_parameters_velu
    init_runtime = algo.field.init_runtime
    kps = algo.formula.kps
    xisog = algo.formula.xisog
    xeval = algo.formula.xeval
    sJ_list = algo.formula.sJ_list
    xmul = algo.curve.xmul
    HYBRID_BOUND = algo.formula.HYBRID_BOUND

    A = [algo.curve.field(2), algo.curve.field(4)]
    if True:

        # T_p belongs to E[pi - 1]
        # T_m belongs to E[pi + 1]
        T_p, T_m = generators(A)

    else:

        # T_m belongs to E[pi - 1]
        # T_p belongs to E[pi + 1]
        T_m, T_p = generators(A)

    assert len(L) == n

    original_stdout = sys.stdout # Save a reference to the original standard output
    multievaluation = {True:'scaled', False:'unscaled'}[setting.multievaluation]
    path = resource_filename('sibc', "data/ijk/" + algo.curve.model + '/csidh/' + algo.curve.name + '-' + multievaluation)
    with open(path, 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        parameters = dict()
        for idx in range(0, n, 1):

            # -------------------------------------------------------------
            # Random kernel point
            Tp = list(T_p)
            for i in range(0, n, 1):
                if i != idx:
                    Tp = xmul(Tp, A, i)

            # -------------------------------------------------------------
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in kps, xisog, and xeval
            if L[idx] == 3:
                b = 0
                c = 0
            else:
                b = int(floor(sqrt(L[idx] - 1) / 2.0))
                c = int(floor((L[idx] - 1.0) / (4.0 * b)))

            b += 1
            parameters[str(idx)] = []
            for j in range(0, int(floor(sqrt(pi * (b - 1)) / 1.0)), 1):

                b -= 1
                c = int(floor((L[idx] - 1.0) / (4.0 * b)))
                set_parameters_velu(b, c, idx)

                total_cost = [0, 0, 0]
                # -------------------------------------------------------------
                # kps procedure
                init_runtime()
                kps(Tp, A, idx)
                
                total_cost[0] += algo.field.fpmul
                total_cost[1] += algo.field.fpsqr
                total_cost[2] += algo.field.fpadd

                # -------------------------------------------------------------
                # xisog
                init_runtime()
                B = xisog(A, idx)
                
                total_cost[0] += algo.field.fpmul
                total_cost[1] += algo.field.fpsqr
                total_cost[2] += algo.field.fpadd

                # -------------------------------------------------------------
                # xeval bench
                init_runtime()
                if L[idx] <= HYBRID_BOUND:
                    Tm = xeval(T_m, idx)
                else:
                    Tm = xeval(T_m, A)
                
                total_cost[0] += algo.field.fpmul
                total_cost[1] += algo.field.fpsqr
                total_cost[2] += algo.field.fpadd

                # assert(validate(B))
                parameters[str(idx)].append(
                    (
                        b,
                        c,
                        [algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd],
                        total_cost[0] + total_cost[1]
                    )
                )

            if L[idx] == 3:
                parameters[str(idx)] = (0, 0, None, None)
            else:
                parameters[str(idx)] = min(
                    parameters[str(idx)], key=lambda tup: tup[3]
                )

            print(parameters[str(idx)][0], parameters[str(idx)][1])
    return attrdict(name='csidh-precompute-parameters', **locals())
