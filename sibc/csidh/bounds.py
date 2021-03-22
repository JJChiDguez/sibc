import click
import numpy

from math import ceil, log

from sibc.common import attrdict, printl
tmp_dir = "./"  # Drop the files in the current working directory

@click.command()
@click.pass_context
def csidh_bounds(ctx):
    """ Greedy-based search of optimal exponents """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    assert setting.uninitialized, 'option -u (--uninitialized) is required!'
    
    L = algo.params.L
    n = algo.params.n
    security = algo.gae.security
    strategy_block_cost = algo.gae.strategy_block_cost
    basis = numpy.eye(n, dtype=int)
    measure = algo.curve.measure

    print_intv = lambda v, n: ', '.join(list(map(format, v, ['2d'] * n)))

    # MITM procedure
    keyspace = 256.00  # For ensuring 128-bits of classical security
    # keyspace = 384.00              # For ensuring 192-bits of classical security
    # vOW Golden Collision Search
    # keyspace = 220.2295746338436   # For ensuring 128-bits of classical security
    # keyspace = 305.5629079671769   # For ensuring 192-bits of classical security

    # The next function computes the set \mu() given in the paper
    def neighboring_intvec(seq_i, L, IN, OUT):
        nonlocal keyspace
        if OUT[2] >= keyspace:
            return OUT

        else:
            minimum = IN
            if measure(IN[1]) >= measure(OUT[1]):

                for j in seq_i:

                    if OUT[0][j] > 0:
                        current_cost, _, _, _, _ = strategy_block_cost(
                            L, OUT[0] + basis[j]
                        )
                        tmp = neighboring_intvec(
                            seq_i,
                            L,
                            IN,
                            (
                                OUT[0] + basis[j],
                                current_cost,
                                security(OUT[0] + basis[j], len(L)),
                            ),
                        )
                        if measure(minimum[1]) >= measure(tmp[1]):
                            minimum = tmp

            return minimum

    # Finally, the next functions is the implementation of algorithm 2.0
    def optimal_bounds(L, b, r):

        assert r >= 1
        n = len(L)

        RNC, _, _, _, _ = strategy_block_cost(L, b)
        SEC = security(b, n)
        e = b

        for i in range(0, n, 1):

            # The algorithm proceed by looking the best bounds when e_i <- e_i - 1
            seq_i = [k for k in range(n) if k != i]
            (e_tmp, RNC_tmp, SEC_tmp) = (e, RNC, SEC)
            if e[i] > r:

                # Set the new possible optimal bounds
                temporal_cost, _, _, _, _ = strategy_block_cost(
                    L, e - r * basis[i]
                )
                (e_tmp, RNC_tmp, SEC_tmp) = neighboring_intvec(
                    seq_i,
                    L,
                    (e, RNC, SEC),
                    (
                        e - r * basis[i],
                        temporal_cost,
                        security(e - r * basis[i], n),
                    ),
                )

            (e, RNC, SEC) = min(
                [(e_tmp, RNC_tmp, SEC_tmp), (e, RNC, SEC)],
                key=lambda t: measure(t[1]),
            )

            print("[Security := %f]" % SEC, end="\t")
            print(
                "decreasing: e_{"
                + print_intv([i], 1)
                + "}"
                + ", and increasing each e_j with j != "
                + print_intv([i], 1)
                + "; current optimal running-time: %7.3f" % measure(RNC)
            )
            print("[" + print_intv(e, n) + "]\n")

        # --------------------------------------------------------------------------------------------------
        f = open(
            tmp_dir
            + "csidh_"
            + setting.prime
            + "_"
            + setting.style
            + "_e"
            + setting.exponent
            + ".py",
            "w",
        )
        #f.write('m = [' + ', '.join([str(ei) for ei in e[::-1]]) + ']')
        f.write(' '.join(map(str, e[::-1])))
        f.close()
        # --------------------------------------------------------------------------------------------------
        return (e, RNC)

    ''' -------------------------------------------------------------------------------------
        Number of degree-(l_i) isogeny constructions to be performed: m_i
        ------------------------------------------------------------------------------------- '''

    # ==========================================================================
    m = int(setting.exponent)

    # ---
    k = 3
    # Next integer vector bount is given in Onuki et al. manuscript
    print(
        "\n_______________________________________________________________________________________________________________________________"
    )
    print("List of small odd primes")
    printl("L", L[::-1], n // k)
    print("\nInitial integer vector of bounts (b_0, ..., b_%d)" % n)

    if m == 1:
        if setting.style == 'wd1' or setting.style == 'df':
            assert n >= int(keyspace), f'not enough prime factors in p + 1: {n}'
            e = [1] * int(keyspace) + [0] * (n - int(keyspace))
        else:
            assert n >= int(ceil(keyspace / log(3,2))), f'not enough prime factors in p + 1: {n}'
            e = [1] * 162 + [0] * (n - 162)

        # --------------------------------------------------------------------------------------------------
        f = open(
            tmp_dir
            + "csidh_"
            + setting.prime
            + "_"
            + setting.style
            + "_e"
            + setting.exponent
            + ".py",
            "w",
        )
        #f.write('m = [' + ', '.join([str(ei) for ei in e]) + ']')
        f.write(' '.join(map(str, e)))
        f.close()
        # --------------------------------------------------------------------------------------------------

    else:
        e = [m] * (n - 1) + [(3 * m) // 2]
        e = numpy.array(e)
        stop = False
        for i in range(0, n, 1):
            for j in range(0, m, 1):
                e = e - basis[i]
                if security(e, n) < keyspace:
                    e = e + basis[i]
                    stop = True
                    break

            if stop:
                break
        printl("e", e, n // k)
        RUNNING_TIME, _, _, _, _ = strategy_block_cost(L[::-1], e)

        print(
            "// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M"
            % (
                RUNNING_TIME[0] / (10.0 ** 6),
                RUNNING_TIME[1] / (10.0 ** 6),
                RUNNING_TIME[2] / (10.0 ** 6),
                measure(RUNNING_TIME) / (10.0 ** 6),
            )
        )
        print("\tSecurity ~ %f\n" % security(e, n))

        print(
            "_______________________________________________________________________________________________________________________________"
        )
        print("We proceed by searching a better integer vector of bounds\n")
        r = 1
        for k in range(1, int(ceil((1.0 * m) / (1.0 * r)))):
            e, RNC = optimal_bounds(L[::-1], e, r)

        print(
            "_______________________________________________________________________________________________________________________________\n"
        )
    return attrdict(name='bounds', **locals())
