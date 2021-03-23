from pkg_resources import resource_filename
import click

from sibc.common import attrdict, geometric_serie, rounds

@click.command()
@click.pass_context
def csidh_main(ctx):
    """ Random instance example of a key-exchange """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
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
    random_exponents = algo.gae.random_exponents
    print_exponents = algo.gae.print_exponents

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
