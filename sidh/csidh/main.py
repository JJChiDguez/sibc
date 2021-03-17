from pkg_resources import resource_filename
import click

from sidh.common import attrdict, geometric_serie, rounds
from sidh.constants import strategy_data

@click.command()
@click.pass_context
def csidh_main(ctx):
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    coeff = algo.curve.coeff
    p = algo.params.p
    L = algo.params.L
    m = algo.params.m
    n = algo.params.n
    SQR, ADD = algo.curve.SQR, algo.curve.ADD
    init_runtime = algo.field.init_runtime
    validate = algo.curve.issupersingular
    measure = algo.curve.measure
    GAE_at_0 = algo.gae.GAE_at_0
    GAE_at_A = algo.gae.GAE_at_A
    strategy_block_cost = algo.gae.strategy_block_cost
    random_exponents = algo.gae.random_exponents
    print_exponents = algo.gae.print_exponents

    if algo.formula.name != 'tvelu':
        set_parameters_velu = algo.formula.set_parameters_velu

    temporal_m = list(set(m))
    if len(temporal_m) > 1:
        # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
        LABEL_m = 'different_bounds'
    else:
        # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
        LABEL_m = 'with_same_bounds'

    if setting.tuned:
        verb = '-suitable'
    else:
        verb = '-classical'

    # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
    m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
    r_out, L_out, R_out = rounds(m_prime[::-1], n)
    for j in range(0, len(r_out), 1):

        R_out[j] = list([L[::-1][k] for k in R_out[j]])
        L_out[j] = list([L[::-1][k] for k in L_out[j]])

    file_path = (
        "data/strategies/"
        + setting.algorithm
        + '-'
        + setting.prime
        + '-'
        + setting.style
        + '-'
        + setting.formula
        + '-'
        + LABEL_m
        + verb
    )
    file_path = resource_filename('sidh', file_path)
    try:
        f = open(file_path)
        print("// Strategies to be read from a file")
        S_out = []
        for i in range(0, len(r_out), 1):

            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            S_out.append(tmp)

        f.close()

    except IOError:

        print("// Strategies to be computed")
        C_out, L_out, R_out, S_out, r_out = strategy_block_cost(
            L[::-1], m[::-1]
        )
        f = open(file_path, 'w')
        for i in range(0, len(r_out)):

            f.writelines(' '.join([str(tmp) for tmp in S_out[i]]) + '\n')

        f.close()

    print(
        "// The running time is assuming S = %1.2f x M and a = %1.2f x M, and giving in millions of field operations.\n"
        % (SQR, ADD)
    )

    ''' -------------------------------------------------------------------------------------
        Main
        ------------------------------------------------------------------------------------- '''
    print_exponents("// Exponent bounds", m)
    print("\n// ===================== \033[0;33mPublic Key Generation\033[0m")

    # ------------------------------------------------------------------------- Alice
    print("// --- \033[0;35mAlice\033[0m")
    init_runtime()
    a_private = random_exponents(m)
    a_public = GAE_at_0(a_private)

    print(
        "// Running time (GAE):\t\t\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            algo.field.fpmul / (10.0 ** 6),
            algo.field.fpsqr / (10.0 ** 6),
            algo.field.fpadd / (10.0 ** 6),
            measure([algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]) / (10.0 ** 6),
        )
    )
    print_exponents("sk_a", a_private)
    print("pk_a := %s;" % coeff(a_public))

    # ------------------------------------------------------------------------- Bob
    print("\n// --- \033[0;34mBob\033[0m")
    init_runtime()
    b_private = random_exponents(m)
    b_public = GAE_at_0(b_private)
    
    print(
        "// Running time (GAE):\t\t\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            algo.field.fpmul / (10.0 ** 6),
            algo.field.fpsqr / (10.0 ** 6),
            algo.field.fpadd / (10.0 ** 6),
            measure([algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]) / (10.0 ** 6),
        )
    )
    print_exponents("sk_b", b_private)
    print("pk_b := %s;" % coeff(b_public))

    print("\n// ===================== \033[0;33mSecret Sharing Computation\033[0m")
    # ------------------------------------------------------------------------- Alice
    print("// --- \033[0;35mAlice\033[0m")
    init_runtime()
    public_validation = validate(b_public)
    assert public_validation
    
    print(
        "// Running time (key validation):\t%2.3fM + %2.3fS + %2.3fa = %2.3fM,"
        % (
            algo.field.fpmul / (10.0 ** 6),
            algo.field.fpsqr / (10.0 ** 6),
            algo.field.fpadd / (10.0 ** 6),
            measure([algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]) / (10.0 ** 6),
        )
    )

    init_runtime()
    ss_a = GAE_at_A(a_private, b_public)
    print(
        "// Running time (GAE + key validation):\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            algo.field.fpmul / (10.0 ** 6),
            algo.field.fpsqr / (10.0 ** 6),
            algo.field.fpadd / (10.0 ** 6),
            measure([algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]) / (10.0 ** 6),
        )
    )
    print("ss_a := %s;\n" % coeff(ss_a))

    # ------------------------------------------------------------------------- Bob
    print("// --- \033[0;34mBob\033[0m")
    init_runtime()
    public_validation = validate(a_public)
    assert public_validation
    
    print(
        "// Running time (key validation):\t%2.3fM + %2.3fS + %2.3fa = %2.3fM,"
        % (
            algo.field.fpmul / (10.0 ** 6),
            algo.field.fpsqr / (10.0 ** 6),
            algo.field.fpadd / (10.0 ** 6),
            measure([algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]) / (10.0 ** 6),
        )
    )

    init_runtime()
    ss_b = GAE_at_A(b_private, a_public)
    
    print(
        "// Running time (GAE + key validation):\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            algo.field.fpmul / (10.0 ** 6),
            algo.field.fpsqr / (10.0 ** 6),
            algo.field.fpadd / (10.0 ** 6),
            measure([algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]) / (10.0 ** 6),
        )
    )
    print("ss_b := %s;" % coeff(ss_b))

    try:
        assert(coeff(ss_a) == coeff(ss_b))
        print('\n\x1b[0;30;43m' + 'Successfully passed!' + '\x1b[0m')
    except:
        raise TypeError(
            '\x1b[0;30;41m'
            + 'Great Scott!... The sky is falling. NOT PASSED!!!'
            + '\x1b[0m'
        )

    return attrdict(name='csidh-main', **locals())
