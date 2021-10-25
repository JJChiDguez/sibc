import click
from math import floor, sqrt
from sibc.common import attrdict
from sibc.math import cswap

from random import SystemRandom

@click.command()
@click.pass_context
def bsidh_test(ctx):
    """ GF(p²)-operation cost of kps, xisog, and xeval """
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

    print("p := 0x%X;" % p)
    print("fp := GF(p);")
    print("_<x> := PolynomialRing(fp);")
    print("fp2<u> := ext<fp | x^2 + 1>;")
    print("Pr<x> := PolynomialRing(fp2);")

    # ---Generators in E[p + 1]
    PA = list([algo.strategy.PA, field(1)])
    QA = list([algo.strategy.QA, field(1)])
    PQA = list([algo.strategy.PQA, field(1)])
    # ---Generators in E[p - 1]
    PB = list([algo.strategy.PB, field(1)])
    QB = list([algo.strategy.QB, field(1)])
    PQB = list([algo.strategy.PQB, field(1)])

    print("E := EllipticCurve(x^3 + (%s) * x^2 + x);" % coeff(A))

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

    print("\n// Verifying torsion-(p + 1) points")
    print("// x([p + 1]PA) = (1:0)?\t", isinfinity(S))
    print("// x([p + 1]QA) = (1:0)?\t", isinfinity(T))
    print("// x([p + 1]PQA) = (1:0)?\t", isinfinity(ST))

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
    print("\n// Verifying torsion-(p - 1) points")
    print("// x([p - 1]PB) = (1:0)?\t", isinfinity(S))
    print("// x([p - 1]QB) = (1:0)?\t", isinfinity(T))
    print("// x([p - 1]PQB) = (1:0)?\t", isinfinity(ST))

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

    print("\n// Verifying orders")
    assert isfullorder(cofactor_multiples(S, A, range(0, np, 1)))
    assert isfullorder(cofactor_multiples(T, A, range(0, np, 1)))
    assert isfullorder(cofactor_multiples(ST, A, range(0, np, 1)))
    print(
        "// PA is a full order point?\t",
        isfullorder(cofactor_multiples(S, A, range(0, np, 1))),
    )
    print(
        "// QA is a full order point?\t",
        isfullorder(cofactor_multiples(T, A, range(0, np, 1))),
    )
    print(
        "// QPA is a full order point?\t",
        isfullorder(cofactor_multiples(ST, A, range(0, np, 1))),
    )

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

    print("\n// Verifying orders")
    assert isfullorder(cofactor_multiples(S, A, range(np, np + nm, 1)))
    assert isfullorder(cofactor_multiples(T, A, range(np, np + nm, 1)))
    assert isfullorder(cofactor_multiples(ST, A, range(np, np + nm, 1)))
    print(
        "// PB is a full order point?\t",
        isfullorder(cofactor_multiples(S, A, range(np, np + nm, 1))),
    )
    print(
        "// QB is a full order point?\t",
        isfullorder(cofactor_multiples(T, A, range(np, np + nm, 1))),
    )
    print(
        "// QPB is a full order point?\t",
        isfullorder(cofactor_multiples(ST, A, range(np, np + nm, 1))),
    )

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
    R = Ladder3pt(k, S, T, ST, algo.curve.field(2))
    T_p = list(R)
    T_m = list(S)
    print(
        "\n// Now, we proceed by performing xisog with input curve equals the output curve of the previous one experiment, and using torsion-(p + 1) points"
    )
    for idx in range(0, np, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(idx + 1, np, 1):
            Tp = xmul(Tp, A, i)

        print("// l:\t%7d |" % L[idx], end="")
        total_cost = [0, 0, 0]

        if setting.formula != 'tvelu':

            if setting.tuned:
                set_parameters_velu(algo.formula.sJ_list[idx], algo.formula.sI_list[idx], idx)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b'
                # from example 4.12 of https://eprint.iacr.org/2020/341 These
                # paramters are required in kps, xisog, and xeval
                if L[idx] <= 4:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(L[idx] - 1) / 2.0))
                    c = int(floor((L[idx] - 1.0) / (4.0 * b)))

                set_parameters_velu(b, c, idx)

            print_parameters_velu()

        # -------------------------------------------------------------
        # kps procedure
        init_runtime_basefield()
        init_runtime_field()
        kps(Tp, A, idx)
        show_runtime_field("kps")
        
        total_cost[0] += algo.field.fp2mul
        total_cost[1] += algo.field.fp2sqr
        total_cost[2] += algo.field.fp2add

        # -------------------------------------------------------------
        # xisog
        init_runtime_basefield()
        init_runtime_field()
        Tp[0], A[0] = cswap(Tp[0], A[0], L[idx] == 4)
        Tp[1], A[1] = cswap(Tp[1], A[1], L[idx] == 4)
        B = xisog(A, idx)
        Tp[0], A[0] = cswap(Tp[0], A[0], L[idx] == 4)
        Tp[1], A[1] = cswap(Tp[1], A[1], L[idx] == 4)
        show_runtime_field("xisog")
        
        total_cost[0] += algo.field.fp2mul
        total_cost[1] += algo.field.fp2sqr
        total_cost[2] += algo.field.fp2add

        # -------------------------------------------------------------
        # xeval: kernel point determined by the next isogeny evaluation
        init_runtime_basefield()
        init_runtime_field()
        if (
            setting.formula == 'tvelu'
            or (setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND)
            or (L[idx] == 4)
        ):
            T_p = xeval(T_p, idx)
        else:
            T_p = xeval(T_p, A)

        # xeval bench
        init_runtime_basefield()
        init_runtime_field()
        if (
            setting.formula == 'tvelu'
            or (setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND)
            or (L[idx] == 4)
        ):
            T_m = xeval(T_m, idx)
        else:
            T_m = xeval(T_m, A)

        show_runtime_field("xeval")
        
        total_cost[0] += algo.field.fp2mul
        total_cost[1] += algo.field.fp2sqr
        total_cost[2] += algo.field.fp2add
        print("|| cost: %8d" % (total_cost[0] + total_cost[1]), end=" ")
        print(
            "|| ratio: %1.3f"
            % ((total_cost[0] + total_cost[1]) / (L[idx] + 2.0))
        )

        # assert(validate(B))
        A = list(B)

        # print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
        # print("assert(Random(B) * (p + 1) eq B!0);")
        # print("BOOL, Q := IsPoint(B, fp!%d/%d);" % (T_m[0], T_m[1]))
        # print("assert(BOOL);")

    print(
        "\n// All the l_i's have been processed, output of xisog corresponds with the given below"
    )
    # print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
    print("B := EllipticCurve(x^3 + (%s) * x^2 + x);" % coeff(A))
    print("assert(Random(B) * (p + 1) eq B!0);")

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
    R = Ladder3pt(k, S, T, ST, algo.curve.field(2))
    T_p = list(R)
    T_m = list(S)
    print(
        "\n// Now, we proceed by performing xisog with input curve equals the output curve of the previous one experiment, and using torsion-(p - 1) points"
    )
    for idx in range(np, np + nm, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(idx + 1, np + nm, 1):
            Tp = xmul(Tp, A, i)

        print("// l:\t%7d |" % L[idx], end="")
        total_cost = [0, 0, 0]

        if setting.formula != 'tvelu':

            if setting.tuned:
                set_parameters_velu(algo.formula.sJ_list[idx], algo.formula.sI_list[idx], idx)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b'
                # from example 4.12 of https://eprint.iacr.org/2020/341 These
                # paramters are required in kps, xisog, and xeval
                if L[idx] == 3:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(L[idx] - 1) / 2.0))
                    c = int(floor((L[idx] - 1.0) / (4.0 * b)))

                set_parameters_velu(b, c, idx)

            print_parameters_velu()

        # -------------------------------------------------------------
        # kps procedure
        init_runtime_basefield()
        init_runtime_field()
        kps(Tp, A, idx)
        show_runtime_field("kps")
        
        total_cost[0] += algo.field.fp2mul
        total_cost[1] += algo.field.fp2sqr
        total_cost[2] += algo.field.fp2add

        # -------------------------------------------------------------
        # xisog
        init_runtime_basefield()
        init_runtime_field()
        B = xisog(A, idx)
        show_runtime_field("xisog")
        
        total_cost[0] += algo.field.fp2mul
        total_cost[1] += algo.field.fp2sqr
        total_cost[2] += algo.field.fp2add

        # -------------------------------------------------------------
        # xeval: kernel point determined by the next isogeny evaluation
        init_runtime_basefield()
        init_runtime_field()
        if setting.formula == 'tvelu' or (
            setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND
        ):
            T_p = xeval(T_p, idx)
        else:
            T_p = xeval(T_p, A)

        # xeval bench
        init_runtime_basefield()
        init_runtime_field()
        if setting.formula == 'tvelu' or (
            setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND
        ):
            T_m = xeval(T_m, idx)
        else:
            T_m = xeval(T_m, A)

        show_runtime_field("xeval")
        
        total_cost[0] += algo.field.fp2mul
        total_cost[1] += algo.field.fp2sqr
        total_cost[2] += algo.field.fp2add
        print("|| cost: %8d" % (total_cost[0] + total_cost[1]), end=" ")
        print(
            "|| ratio: %1.3f"
            % ((total_cost[0] + total_cost[1]) / (L[idx] + 2.0))
        )

        # assert(validate(B))
        A = list(B)

        # print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
        # print("assert(Random(B) * (p + 1) eq B!0);")
        # print("BOOL, Q := IsPoint(B, fp!%d/%d);" % (T_m[0], T_m[1]))
        # print("assert(BOOL);")

    print(
        "\n// All the l_i's have been processed, output of xisog corresponds with the given below"
    )
    
    print("B := EllipticCurve(x^3 + (%s) * x^2 + x);" % coeff(A))
    print("assert(Random(B) * (p + 1) eq B!0);")
    # """

    print(
        "\n\"If no errors were showed using magma calculator, then all experiments were successfully passed!\";"
    )
    print("// copy and paste it at http://magma.maths.usyd.edu.au/calc/\n")
    return attrdict(name='bsidh-test', **locals())
