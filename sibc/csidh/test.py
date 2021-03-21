import click
from math import floor, sqrt

from sibc.common import attrdict

@click.command()
@click.pass_context
def csidh_test(ctx):
    """ GF(p)-operation cost of kps, xisog, and xeval """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    generators = algo.curve.generators
    coeff = algo.curve.coeff
    p = algo.params.p
    L = algo.params.L
    n = algo.params.n
    A = [ algo.curve.field(2), algo.curve.field(4)]
    xmul = algo.curve.xmul
    if algo.formula.name != 'tvelu':
        set_parameters_velu = algo.formula.set_parameters_velu
        print_parameters_velu = algo.formula.print_parameters_velu
        HYBRID_BOUND = algo.formula.HYBRID_BOUND
    
    init_runtime = algo.field.init_runtime
    kps = algo.formula.kps
    show_runtime = algo.field.show_runtime
    xisog = algo.formula.xisog
    xeval = algo.formula.xeval

    total_cost = [0, 0, 0]
    print("p := 0x%X;" % p)
    print("fp := GF(p);")
    print("P<x> := PolynomialRing(fp);")

    print("E := EllipticCurve(x^3 + (%s) * x^2 + x);" % coeff(A))

    if True:

        # T_p belongs to E[pi - 1]
        # T_m belongs to E[pi + 1]
        T_p, T_m = generators(A)

    else:

        # T_m belongs to E[pi - 1]
        # T_p belongs to E[pi + 1]
        T_m, T_p = generators(A)

    assert len(L) == n
    print(
        "// Now, we proceed by performing xisog with input curve equals the output curve of the previous one experiment."
    )
    for idx in range(0, n, 1):

        # -------------------------------------------------------------
        # Random kernel point
        Tp = list(T_p)
        for i in range(idx + 1, n, 1):
            Tp = xmul(Tp, A, i)

        print("// l:\t%7d |" % L[idx], end="")
        total_cost = [0, 0, 0]

        if setting.formula != 'tvelu':

            if setting.tuned:
                set_parameters_velu(algo.formula.sJ_list[idx], algo.formula.sI_list[idx], idx)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in kps, xisog, and xeval
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
        init_runtime()
        kps(Tp, A, idx)
        show_runtime("kps")
        total_cost[0] += algo.field.fpmul
        total_cost[1] += algo.field.fpsqr
        total_cost[2] += algo.field.fpadd

        # -------------------------------------------------------------
        # xisog
        init_runtime()
        B = xisog(A, idx)
        show_runtime("xisog")
        total_cost[0] += algo.field.fpmul
        total_cost[1] += algo.field.fpsqr
        total_cost[2] += algo.field.fpadd

        # -------------------------------------------------------------
        # xeval: kernel point determined by the next isogeny evaluation
        init_runtime()
        if setting.formula == 'tvelu' or (
            setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND
        ):
            T_p = xeval(T_p, idx)
        else:
            T_p = xeval(T_p, A)

        # xeval bench
        init_runtime()
        if setting.formula == 'tvelu' or (
            setting.formula == 'hvelu' and L[idx] <= HYBRID_BOUND
        ):
            T_m = xeval(T_m, idx)
        else:
            T_m = xeval(T_m, A)

        show_runtime("xeval")
        total_cost[0] += algo.field.fpmul
        total_cost[1] += algo.field.fpsqr
        total_cost[2] += algo.field.fpadd
        print("|| cost: %7d" % (total_cost[0] + total_cost[1]), end=" ")
        print(
            "|| ratio: %1.3f"
            % ((total_cost[0] + total_cost[1]) / (L[idx] + 2.0))
        )

        # assert(validate(B))
        A = list(B)

        #print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
        #print("assert(Random(B) * (p + 1) eq B!0);")
        #print("BOOL, Q := IsPoint(B, fp!%d/%d);" % (T_m[0], T_m[1]))
        #print("assert(BOOL);")

    print(
        "\n// All the l_i's have been processed, output of xisog corresponds with the given below"
    )
    print("B := EllipticCurve(x^3 + (%s) * x^2 + x);" % coeff(A))
    print("assert(Random(B) * (p + 1) eq B!0);")

    print(
        "\n\"If no errors were showed using magma calculator, then all experiments were successfully passed!\";"
    )
    print("// copy and paste it at http://magma.maths.usyd.edu.au/calc/\n")
    return attrdict(name='csidh-test', **locals())
