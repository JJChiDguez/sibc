import sys
import click
from pkg_resources import resource_filename

from math import floor, sqrt, pi
from sibc.common import attrdict
from sibc.math import cswap

from random import SystemRandom

@click.command()
@click.pass_context
def bsidh_precompute_parameters(ctx):
    """ Precomputation of tuned velusqrt parameters """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    p = algo.params.p
    np = algo.params.np
    Ep = algo.params.Ep
    nm = algo.params.nm
    Em = algo.params.Em
    Lp = algo.params.Lp
    Lm = algo.params.Lm
    L = list(Lp + Lm)

    A = [ algo.curve.field(8), algo.curve.field(4)]
    xmul = algo.curve.xmul
    isinfinity = algo.curve.isinfinity
    coeff = algo.curve.coeff
    isfullorder = algo.curve.isfullorder
    cofactor_multiples = algo.curve.cofactor_multiples
    Ladder3pt = algo.curve.Ladder3pt
    if algo.formula.name != 'tvelu':
        set_parameters_velu = algo.formula.set_parameters_velu
        print_parameters_velu = algo.formula.print_parameters_velu
        HYBRID_BOUND = algo.formula.HYBRID_BOUND
    
    init_runtime_field = algo.field.init_runtime
    show_runtime_field = algo.field.show_runtime
    init_runtime_basefield = algo.basefield.init_runtime
    show_runtime_basefield = algo.basefield.show_runtime
    kps = algo.formula.kps
    xisog = algo.formula.xisog
    xeval = algo.formula.xeval

    field = algo.field
    random = SystemRandom()

    # ---Generators in E[p + 1]
    PA = [algo.strategy.PA, field(1)]
    QA = [algo.strategy.QA, field(1)]
    PQA = [algo.strategy.PQA, field(1)]
    # ---Generators in E[p - 1]
    PB = [algo.strategy.PB, field(1)]
    QB = [algo.strategy.QB, field(1)]
    PQB = [algo.strategy.PQB, field(1)]

    S = list(PA)
    T = list(QA)
    ST = list(PQA)

    for i in range(0, np, 1):
        for idx in range(0, Ep[i], 1):
            S = xmul(S, A, i)
            T = xmul(T, A, i)
            ST = xmul(ST, A, i)

    assert isinfinity(S)
    assert isinfinity(T)
    assert isinfinity(ST)

    S = list(PB)
    T = list(QB)
    ST = list(PQB)

    for i in range(np, np + nm, 1):
        for idx in range(0, Em[i - np], 1):
            S = xmul(S, A, i)
            T = xmul(T, A, i)
            ST = xmul(ST, A, i)

    assert isinfinity(S)
    assert isinfinity(T)
    assert isinfinity(ST)

    # Case (p + 1)
    S = list(PA)
    T = list(QA)
    ST = list(PQA)

    assert isinfinity(S) == False
    assert isinfinity(T) == False
    assert isinfinity(ST) == False

    for i in range(0, np, 1):
        for idx in range(0, Ep[i] - 1, 1):
            S = xmul(S, A, i)
            T = xmul(T, A, i)
            ST = xmul(ST, A, i)

    assert isfullorder(cofactor_multiples(S, A, range(0, np, 1)))
    assert isfullorder(cofactor_multiples(T, A, range(0, np, 1)))
    assert isfullorder(cofactor_multiples(ST, A, range(0, np, 1)))

    # Case (p - 1)
    S = list(PB)
    T = list(QB)
    ST = list(PQB)

    assert isinfinity(S) == False
    assert isinfinity(T) == False
    assert isinfinity(ST) == False

    for i in range(np, np + nm, 1):
        for idx in range(0, Em[i - np] - 1, 1):
            S = xmul(S, A, i)
            T = xmul(T, A, i)
            ST = xmul(ST, A, i)

    assert isfullorder(cofactor_multiples(S, A, range(np, np + nm, 1)))
    assert isfullorder(cofactor_multiples(T, A, range(np, np + nm, 1)))
    assert isfullorder(cofactor_multiples(ST, A, range(np, np + nm, 1)))

    # Three point ladder: case (p + 1)
    S = list(PA)
    T = list(QA)
    ST = list(PQA)

    assert isinfinity(S) == False
    assert isinfinity(T) == False

    for i in range(0, np, 1):
        for idx in range(0, Ep[i] - 1, 1):
            S = xmul(S, A, i)
            T = xmul(T, A, i)
            ST = xmul(ST, A, i)

    k = random.randint(0, p)
    R = Ladder3pt(k, S, T, ST, algo.curve.field(6))
    T_p = list(R)
    T_m = list(S)

    original_stdout = sys.stdout # Save a reference to the original standard output
    multievaluation = {True:'scaled', False:'unscaled'}[setting.multievaluation]
    path = resource_filename('sibc', "data/ijk/" + algo.curve.model + '/bsidh/' + algo.curve.name + '-' + multievaluation)
    with open(path, 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        parameters = dict()
        for idx in range(0, np, 1):

            # -------------------------------------------------------------
            # Random kernel point
            Tp = list(T_p)
            for i in range(0, np, 1):
                if i != idx:
                    Tp = xmul(Tp, A, i)

            # -------------------------------------------------------------
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in kps, xisog, and xEVAL
            if L[idx] <= 4:
                b = 0
                c = 0
            else:
                b = int(floor(sqrt(L[idx] - 1) / 2.0))
                c = int(floor((L[idx] - 1.0) / (4.0 * b)))

            parameters[str(idx)] = []
            b += 1
            for j in range(0, int(floor(sqrt(pi * (b - 1)) / 1.0)), 1):

                b -= 1
                c = int(floor((L[idx] - 1.0) / (4.0 * b)))
                set_parameters_velu(b, c, idx)

                total_cost = [0, 0, 0]

                # -------------------------------------------------------------
                # kps procedure
                init_runtime_basefield()
                kps(Tp, A, idx)
                
                total_cost[0] += algo.basefield.fpmul
                total_cost[1] += algo.basefield.fpsqr
                total_cost[2] += algo.basefield.fpadd

                # -------------------------------------------------------------
                # xisog
                init_runtime_basefield()
                Tp[0], A[0] = cswap(Tp[0], A[0], L[idx] == 4)
                Tp[1], A[1] = cswap(Tp[1], A[1], L[idx] == 4)
                B = xisog(A, idx)
                Tp[0], A[0] = cswap(Tp[0], A[0], L[idx] == 4)
                Tp[1], A[1] = cswap(Tp[1], A[1], L[idx] == 4)
                
                total_cost[0] += algo.basefield.fpmul
                total_cost[1] += algo.basefield.fpsqr
                total_cost[2] += algo.basefield.fpadd

                # -------------------------------------------------------------
                # xEVAL bench
                init_runtime_basefield()
                if setting.formula == 'tvelu' or (
                    setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND
                    ):
                    Tm = xeval(T_m, idx)
                else:
                    Tm = xeval(T_m, A)
                
                total_cost[0] += algo.basefield.fpmul
                total_cost[1] += algo.basefield.fpsqr
                total_cost[2] += algo.basefield.fpadd

                parameters[str(idx)].append((
                    b,
                    c,
                    [algo.basefield.fpmul, algo.basefield.fpsqr, algo.basefield.fpadd],
                    total_cost[0] + total_cost[1]
                    )
                )

            if L[idx] <= 4:
                parameters[str(idx)] = (0, 0, None, None)
            else:
                parameters[str(idx)] = min(
                    parameters[str(idx)], key=lambda tup: tup[3]
                )

            print(parameters[str(idx)][0], parameters[str(idx)][1])

        # --------------------------------------------------------------------------------------------------------------------
        A = [ algo.curve.field(8), algo.curve.field(4)]
        # Three point ladder: case (p - 1)
        S = list(PB)
        T = list(QB)
        ST = list(PQB)

        assert isinfinity(S) == False
        assert isinfinity(T) == False

        for i in range(np, np + nm, 1):
            for idx in range(0, Em[i - np] - 1, 1):
                S = xmul(S, A, i)
                T = xmul(T, A, i)
                ST = xmul(ST, A, i)

        k = random.randint(0, p)
        R = Ladder3pt(k, S, T, ST, algo.curve.field(6))
        T_p = list(R)
        T_m = list(S)
        for idx in range(np, np + nm, 1):

            # -------------------------------------------------------------
            # Random kernel point
            Tp = list(T_p)
            for i in range(np, np + nm, 1):
                if i != idx:
                    Tp = xmul(Tp, A, i)

            # -------------------------------------------------------------
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in kps, xisog, and xEVAL
            if L[idx] == 3:
                b = 0
                c = 0
            else:
                b = int(floor(sqrt(L[idx] - 1) / 2.0))
                c = int(floor((L[idx] - 1.0) / (4.0 * b)))

            parameters[str(idx)] = []
            b += 1
            for j in range(0, int(floor(sqrt(pi * (b - 1)) / 1.0)), 1):

                b -= 1
                c = int(floor((L[idx] - 1.0) / (4.0 * b)))
                set_parameters_velu(b, c, idx)

                total_cost = [0, 0, 0]
                # -------------------------------------------------------------
                # kps procedure
                init_runtime_basefield()
                kps(Tp, A, idx)
                
                total_cost[0] += algo.basefield.fpmul
                total_cost[1] += algo.basefield.fpsqr
                total_cost[2] += algo.basefield.fpadd

                # -------------------------------------------------------------
                # xisog
                init_runtime_basefield()
                B = xisog(A, idx)
                
                total_cost[0] += algo.basefield.fpmul
                total_cost[1] += algo.basefield.fpsqr
                total_cost[2] += algo.basefield.fpadd

                # -------------------------------------------------------------
                # xEVAL bench
                init_runtime_basefield()
                if setting.formula == 'tvelu' or (
                    setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND
                ):
                    Tm = xeval(T_m, idx)
                else:
                    Tm = xeval(T_m, A)
                
                total_cost[0] += algo.basefield.fpmul
                total_cost[1] += algo.basefield.fpsqr
                total_cost[2] += algo.basefield.fpadd

                parameters[str(idx)].append((
                    b,
                    c,
                    [algo.basefield.fpmul, algo.basefield.fpsqr, algo.basefield.fpadd],
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
        sys.stdout = original_stdout # Reset the standard output to its original value

    return attrdict(name='bsidh-precompute-parameters', **locals())
