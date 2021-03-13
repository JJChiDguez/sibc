import click
from sympy import symbols, floor, sqrt, sign
from pkg_resources import resource_filename
from random import SystemRandom
from math import pi

from sidh.common import attrdict

@click.command()
@click.pass_context
def bsidh_parameters(ctx):
    # This is only valid for svelu and hvelu - if we are called from tvelu, we
    # should exit!
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    tuned_name = ('-classical','-suitable')[setting.tuned]
    curve = algo.curve
    formula = algo.formula
    coeff = curve.coeff
    random = SystemRandom()
    p = algo.params.p
    np = algo.params.np
    Ep = algo.params.Ep
    nm = algo.params.nm
    Em = algo.params.Em
    global_L = algo.curve.L



    # Reading public generators points
    f = open(resource_filename('sidh', 'data/gen/'+ algo.prime))

    # x(PA), x(QA) and x(PA - QA)
    PQA = f.readline()
    PQA = [int(x, 16) for x in PQA.split()]
    PA = [list(PQA[0:2]), [0x1, 0x0]]
    QA = [list(PQA[2:4]), [0x1, 0x0]]
    PQA = [list(PQA[4:6]), [0x1, 0x0]]

    # x(PB), x(QB) and x(PB - QB)
    PQB = f.readline()
    PQB = [int(x, 16) for x in PQB.split()]
    PB = [list(PQB[0:2]), [0x1, 0x0]]
    QB = [list(PQB[2:4]), [0x1, 0x0]]
    PQB = [list(PQB[4:6]), [0x1, 0x0]]

    f.close()

    A = [[0x8, 0x0], [0x4, 0x0]]
    a = coeff(A)

    S = [list(PA[0]), list(PA[1])]
    T = [list(QA[0]), list(QA[1])]
    ST = [list(PQA[0]), list(PQA[1])]

    for i in range(0, np, 1):
        for idx in range(0, Ep[i], 1):
            S = curve.xMUL(S, A, i)
            T = curve.xMUL(T, A, i)
            ST = curve.xMUL(ST, A, i)

    assert curve.isinfinity(S)
    assert curve.isinfinity(T)
    assert curve.isinfinity(ST)

    S = [list(PB[0]), list(PB[1])]
    T = [list(QB[0]), list(QB[1])]
    ST = [list(PQB[0]), list(PQB[1])]

    for i in range(np, np + nm, 1):
        for idx in range(0, Em[i - np], 1):
            S = curve.xMUL(S, A, i)
            T = curve.xMUL(T, A, i)
            ST = curve.xMUL(ST, A, i)

    assert curve.isinfinity(S)
    assert curve.isinfinity(T)
    assert curve.isinfinity(ST)

    # Case (p + 1)
    S = [list(PA[0]), list(PA[1])]
    T = [list(QA[0]), list(QA[1])]
    ST = [list(PQA[0]), list(PQA[1])]

    assert curve.isinfinity(S) == False
    assert curve.isinfinity(T) == False
    assert curve.isinfinity(ST) == False

    for i in range(0, np, 1):
        for idx in range(0, Ep[i] - 1, 1):
            S = curve.xMUL(S, A, i)
            T = curve.xMUL(T, A, i)
            ST = curve.xMUL(ST, A, i)

    # Case (p - 1)
    S = [list(PB[0]), list(PB[1])]
    T = [list(QB[0]), list(QB[1])]
    ST = [list(PQB[0]), list(PQB[1])]

    assert curve.isinfinity(S) == False
    assert curve.isinfinity(T) == False
    assert curve.isinfinity(ST) == False

    for i in range(np, np + nm, 1):
        for idx in range(0, Em[i - np] - 1, 1):
            S = curve.xMUL(S, A, i)
            T = curve.xMUL(T, A, i)
            ST = curve.xMUL(ST, A, i)

    # --------------------------------------------------------------------------------------------------------------------
    # Three point ladder: case (p + 1)
    S = [list(PA[0]), list(PA[1])]
    T = [list(QA[0]), list(QA[1])]
    ST = [list(PQA[0]), list(PQA[1])]

    assert curve.isinfinity(S) == False
    assert curve.isinfinity(T) == False

    for i in range(0, np, 1):
        for idx in range(0, Ep[i] - 1, 1):
            S = curve.xMUL(S, A, i)
            T = curve.xMUL(T, A, i)
            ST = curve.xMUL(ST, A, i)

    k = random.randint(0, p)
    R = curve.Ladder3pt(k, S, T, ST, A)
    # print("k := 0x%X;" % k)
    # print("boolR, R := IsPoint(E, (0x%X + i * 0x%X) / (0x%X + i * 0x%X));" % (R[0][0], R[0][1], R[1][0], R[1][1]))

    T_p = [list(R[0]), list(R[1])]
    T_m = [list(S[0]), list(S[1])]

    parameters = dict()
    for idx in range(0, np, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(0, np, 1):
            if i != idx:
                Tp = curve.xMUL(Tp, A, i)

        # -------------------------------------------------------------
        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
        # These paramters are required in formula.KPs, formula.xISOG, and xEVAL
        if global_L[idx] <= 4:
            b = 0
            c = 0
        else:
            b = int(floor(sqrt(global_L[idx] - 1) / 2.0))
            c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))

        parameters[str(idx)] = []
        b += 1
        for j in range(0, int(floor(sqrt(pi * (b - 1)) / 1.0)), 1):

            b -= 1
            c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))
            formula.set_parameters_velu(b, c, idx)

            total_cost = [0, 0, 0]

            # -------------------------------------------------------------
            # formula.KPs procedure
            formula.fp.fp.set_zero_ops()
            formula.KPs(Tp, A, idx)
            t = formula.fp.fp.get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # -------------------------------------------------------------
            # formula.xISOG
            formula.fp.fp.set_zero_ops()
            Tp[0], A[0] = formula.fp.fp2_cswap(Tp[0], A[0], global_L[idx] == 4)
            Tp[1], A[1] = formula.fp.fp2_cswap(Tp[1], A[1], global_L[idx] == 4)
            B = formula.xISOG(A, idx)
            Tp[0], A[0] = formula.fp.fp2_cswap(Tp[0], A[0], global_L[idx] == 4)
            Tp[1], A[1] = formula.fp.fp2_cswap(Tp[1], A[1], global_L[idx] == 4)
            t = formula.fp.fp.get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # -------------------------------------------------------------
            # xEVAL bench
            formula.fp.fp.set_zero_ops()
            Tm = formula.xEVAL(T_m, A)
            t = formula.fp.fp.get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            parameters[str(idx)].append((b, c, t, total_cost[0] + total_cost[1]))

        if global_L[idx] <= 4:
            parameters[str(idx)] = (0, 0, None, None)
        else:
            parameters[str(idx)] = min(
                parameters[str(idx)], key=lambda tup: tup[3]
            )

        print(parameters[str(idx)][0], parameters[str(idx)][1])

    # --------------------------------------------------------------------------------------------------------------------
    A = [[0x8, 0x0], [0x4, 0x0]]
    # Three point ladder: case (p - 1)
    S = [list(PB[0]), list(PB[1])]
    T = [list(QB[0]), list(QB[1])]
    ST = [list(PQB[0]), list(PQB[1])]

    assert curve.isinfinity(S) == False
    assert curve.isinfinity(T) == False

    for i in range(np, np + nm, 1):
        for idx in range(0, Em[i - np] - 1, 1):
            S = curve.xMUL(S, A, i)
            T = curve.xMUL(T, A, i)
            ST = curve.xMUL(ST, A, i)

    k = random.randint(0, p)
    R = curve.Ladder3pt(k, S, T, ST, A)
    # print("k := 0x%X;" % k)
    # print("boolR, R := IsPoint(E, (0x%X + i * 0x%X) / (0x%X + i * 0x%X));" % (R[0][0], R[0][1], R[1][0], R[1][1]))

    T_p = [list(R[0]), list(R[1])]
    T_m = [list(S[0]), list(S[1])]
    for idx in range(np, np + nm, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(np, np + nm, 1):
            if i != idx:
                Tp = curve.xMUL(Tp, A, i)

        # -------------------------------------------------------------
        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
        # These paramters are required in formula.KPs, formula.xISOG, and xEVAL
        if global_L[idx] == 3:
            b = 0
            c = 0
        else:
            b = int(floor(sqrt(global_L[idx] - 1) / 2.0))
            c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))

        parameters[str(idx)] = []
        b += 1
        for j in range(0, int(floor(sqrt(pi * (b - 1)) / 1.0)), 1):

            b -= 1
            c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))
            formula.set_parameters_velu(b, c, idx)

            total_cost = [0, 0, 0]
            # -------------------------------------------------------------
            # formula.KPs procedure
            formula.fp.fp.set_zero_ops()
            formula.KPs(Tp, A, idx)
            t = formula.fp.fp.get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # -------------------------------------------------------------
            # formula.xISOG
            formula.fp.fp.set_zero_ops()
            B = formula.xISOG(A, idx)
            t = formula.fp.fp.get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            # -------------------------------------------------------------
            # xEVAL bench
            formula.fp.fp.set_zero_ops()
            Tm = formula.xEVAL(T_m, A)
            t = formula.fp.fp.get_ops()
            total_cost[0] += t[0]
            total_cost[1] += t[1]
            total_cost[2] += t[2]

            parameters[str(idx)].append((b, c, t, total_cost[0] + total_cost[1]))

        if global_L[idx] == 3:
            parameters[str(idx)] = (0, 0, None, None)
        else:
            parameters[str(idx)] = min(
                parameters[str(idx)], key=lambda tup: tup[3]
            )

        print(parameters[str(idx)][0], parameters[str(idx)][1])
