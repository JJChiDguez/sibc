from framework import *

if setting.style == 'wd1':
    # Using one torsion point and dummy isogeny constructions
    from csidh.gae_wd1 import *
    delta = 1
elif setting.style == 'wd2':
    # Using two torsion points and dummy isogeny constructions
    from csidh.gae_wd2 import *
    delta = 2
elif setting.style == 'df':
    # Dummy-free approach by using two torsion points
    from csidh.gae_df import *
    delta = 1
else:

    print("  ,-~~-.___.          ")
    print(" / |  '     \\          SYNTAX ERROR ..., run python3 main.py -h for help") 
    print("(  )         0        ")  
    print(" \_/-, ,----'         ")         
    print("    ====           // ")
    print("   /  \-'~;    /~~~(O)")
    print("  /  __/~|   /       |")   
    print("=(  _____| (_________|")
    exit(7)

# -----------------------------------------------------------------
print_intv = lambda v, n: ', '.join(list(map(format, v, ['2d']*n)))

# The next function computes the set \mu() given in the paper
def neighboring_intvec(seq_i, L, IN, OUT):
    if OUT[2] >= 256.00:
        return OUT
    
    else:
        minimum = IN
        if measure(IN[1]) >= measure(OUT[1]):

            for j in seq_i:
            
                    current_cost, _, _, _, _ = strategy_block_cost(L, OUT[0] + basis[j])
                    tmp = neighboring_intvec(   seq_i, 
                                                L, 
                                                IN, 
                                                (OUT[0] + basis[j], current_cost, security(OUT[0] + basis[j], len(L)))
                                            )
                    if measure(minimum[1]) >= measure(tmp[1]):
                        minimum = tmp
        
        return minimum

# Finally, the next functions is the implementation of algorithm 2.0
def optimal_bounds(L, b, r):

    assert(r >= 1)
    n = len(L)

    RNC, _, _, _, _ = strategy_block_cost(L, b)
    SEC = security(b, n)
    e = b

    #for i in range(n-1,-1,-1):
    for i in range(n):

        # The algorithm proceed by looking the best bounds when e_i <- e_i - 1
        seq_i = [ k for k in range(n) if k != i ]
        (e_tmp, RNC_tmp, SEC_tmp) = (e, RNC, SEC)
        if e[i] > r:
        
            # Set the new possible optimal bounds
            temporal_cost, _, _, _, _ = strategy_block_cost(L, e - r*basis[i])
            (e_tmp, RNC_tmp, SEC_tmp) = neighboring_intvec( seq_i,
                                                            L, 
                                                            (e, RNC, SEC),
                                                            (e - r*basis[i], temporal_cost, security(e - r*basis[i], n))
                                        )

        (e, RNC, SEC) = min([(e_tmp, RNC_tmp, SEC_tmp), (e, RNC, SEC)], key=lambda t: measure(t[1]) )
        # --------------------------------------------------------------------------------------------------
        f = open("./tmp/" + setting.style + ".out", "a")
        f.write("decreasing: e_{" + print_intv([i], 1) + "}" +\
                ", and increasing each e_j with j != " + print_intv([i], 1) + "; current optimal running-time: %7.3f\n" % measure(RNC))
        f.write("[" + print_intv(e, n) + "]\n")
        f.close()
        # --------------------------------------------------------------------------------------------------

        print("decreasing: e_{" + print_intv([i], 1) + "}" +\
                ", and increasing each e_j with j != " + print_intv([i], 1) + "; current optimal running-time: %7.3f" % measure(RNC))
        print("[" + print_intv(e, n) + "]\n")
    
    print("Security in bits ~ %f\n" % SEC)
    return (e, RNC)
    
''' -------------------------------------------------------------------------------------
    Number of degree-(l_i) isogeny constructions to be performed: m_i
    ------------------------------------------------------------------------------------- '''

# ==========================================================================
m = 5
 
# ---
k = 3
# Next integer vector bount is given in Onuki et al. manuscript
print("\n_______________________________________________________________________________________________________________________________")
print("List of small odd primes")
printl("L", L[::-1], n // k)
print("\nInitial integer vector of bounts (b_0, ..., b_%d)" % n)

e = [m] * (n - 1) + [(3 * m) // 2]
e = numpy.array(e)
stop = False
for i in range(0, n, 1):

    for j in range(0, m, 1):
        e = e - basis[i]
        if security(e, n) < 256.0:
            e = e + basis[i]
            stop = True
            break

    if stop:
        break
printl("e", e, n // k)
RUNNING_TIME, _, _, _, _ = strategy_block_cost(L[::-1], e)

print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
print("\tSecurity ~ %f\n" % security(e, n))

print("_______________________________________________________________________________________________________________________________")
print("We proceed by searching a better integer vector of bounds\n")
f = open("./tmp/" + setting.style + ".out", "w")
f.write("\n")
f.close()

r = 1
for k in range(1, int(ceil( (1.0*m) / (1.0*r) ))):
    
    e, RNC = optimal_bounds(L[::-1], e, r)
    f = open("./tmp/" + setting.style + ".out", "a")
    f.write("\n")
    f.close()

print("_______________________________________________________________________________________________________________________________\n")
