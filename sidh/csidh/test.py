import click
from sympy import symbols, floor, sqrt, sign

from sidh.common import attrdict

@click.command()
@click.pass_context
def csidh_test(ctx):
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    full_torsion_points = algo.curve.full_torsion_points
    coeff = algo.curve.coeff
    p = algo.params.p
    global_L = L = algo.params.L
    n = algo.params.n
    A = algo.curve.A
    xMUL = algo.curve.xMUL
    set_parameters_velu = algo.formula.set_parameters_velu
    print_parameters_velu = algo.formula.print_parameters_velu
    get_ops = algo.fp.get_ops
    set_zero_ops = algo.fp.set_zero_ops
    sI = algo.formula.sI
    KPs = algo.formula.KPs
    show_ops = algo.fp.show_ops
    xISOG = algo.formula.xISOG
    xEVAL = algo.formula.xEVAL
    HYBRID_BOUND = algo.formula.HYBRID_BOUND

    total_cost = [0, 0, 0]
    print("p := 0x%X;" % p)
    print("fp := GF(p);")
    print("P<x> := PolynomialRing(fp);")

    print("E := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))

    if True:

        # T_p belongs to E[pi - 1]
        # T_m belongs to E[pi + 1]
        T_p, T_m = full_torsion_points(A)

    else:

        # T_m belongs to E[pi - 1]
        # T_p belongs to E[pi + 1]
        T_m, T_p = full_torsion_points(A)

    assert len(L) == n
    print(
        "// Now, we proceed by performing xISOG with input curve equals the output curve of the previous one experiment."
    )
    for idx in range(0, n, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(idx + 1, n, 1):
            Tp = xMUL(Tp, A, i)

        print("// l:\t%7d |" % global_L[idx], end="")
        total_cost = [0, 0, 0]

        if setting.formula != 'tvelu':

            if setting.verbose:
                set_parameters_velu(sJ_list[idx], sI_list[idx], idx)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in KPs, xISOG, and xEVAL
                if global_L[idx] == 3:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(global_L[idx] - 1) / 2.0))
                    c = int(floor((global_L[idx] - 1.0) / (4.0 * b)))

                set_parameters_velu(b, c, idx)

            print_parameters_velu()

        # -------------------------------------------------------------
        # KPs procedure
        set_zero_ops()
        KPs(Tp, A, idx)
        show_ops("Kps", 1.0, 0.0, False)
        t = get_ops()
        total_cost[0] += t[0]
        total_cost[1] += t[1]
        total_cost[2] += t[2]

        # -------------------------------------------------------------
        # xISOG
        set_zero_ops()
        B = xISOG(A, idx)
        show_ops("xISOG", 1.0, 0.0, False)
        t = get_ops()
        total_cost[0] += t[0]
        total_cost[1] += t[1]
        total_cost[2] += t[2]

        # -------------------------------------------------------------
        # xEVAL: kernel point determined by the next isogeny evaluation
        set_zero_ops()
        if setting.formula == 'tvelu' or (
            setting.formula == 'hvelu' and global_L[idx] <= HYBRID_BOUND
        ):
            T_p = xEVAL(T_p, idx)
        else:
            T_p = xEVAL(T_p, A)

        # xEVAL bench
        set_zero_ops()
        if setting.formula == 'tvelu' or (
            setting.formula == 'hvelu' and global_L[idx] <= HYBRID_BOUND
        ):
            T_m = xEVAL(T_m, idx)
        else:
            T_m = xEVAL(T_m, A)

        show_ops("xEVAL", 1.0, 0.0, False)
        t = get_ops()
        total_cost[0] += t[0]
        total_cost[1] += t[1]
        total_cost[2] += t[2]
        print("|| cost: %7d" % (total_cost[0] + total_cost[1]), end=" ")
        print(
            "|| ratio: %1.3f"
            % ((total_cost[0] + total_cost[1]) / (global_L[idx] + 2.0))
        )

        # assert(validate(B))
        A = list(B)

        # print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
        # print("assert(Random(B) * (p + 1) eq B!0);")
        # print("BOOL, Q := IsPoint(B, fp!%d/%d);" % (T_m[0], T_m[1]))
        # print("assert(BOOL);")

    print(
        "\n// All the l_i's have been processed, output of xISOG corresponds with the given below"
    )
    print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
    print("assert(Random(B) * (p + 1) eq B!0);")

    print(
        "\n\"If no errors were showed using magma calculator, then all experiments were successful passed!\";"
    )
    print("// copy and paste it at http://magma.maths.usyd.edu.au/calc/\n")
    return attrdict(name='csidh-test', **locals())

