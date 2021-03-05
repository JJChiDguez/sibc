from sidh.csidh.poly_mul import *

# This function computes the reciprocal of a given polynomial
def reciprocal(f, flen, n):

    if flen < n:
        # truncated reciprocal
        return reciprocal(f + [ 0 ]*(n - flen), n, n)

    if n == 0:
        # Super special base case
        return [1], 1

    if n == 1:

        # Base case when (1/f) mod x^n = 1
        return [1], f[0]

    elif n == 2:

        # Base case when (1/f) = f_0 - f_1*x
        return  [f[0], fp_sub(0, f[1])], fp_sqr(f[0])
    
    elif n == 3:
        
        g, a = reciprocal(f[:2], 2, 2)
        
        t0 = fp_sqr(f[1])
        t1 = fp_mul(f[0], f[2])
        t2 = fp_sub(t1, t0)
        t2 = fp_mul(t2, f[0])
        return [ fp_mul(g[0], a), fp_mul(g[1], a), fp_sub(0, t2) ], fp_sqr(a)
         
    elif n == 4:

        # This case gives the same cost as the general case
        g, a = reciprocal(f[:2], 2, 2)
        t0 = fp_sqr(f[1])
        t1 = fp_mul(g[0], f[2])
        t2 = fp_mul(g[0], f[3])
        t3 = fp_mul(g[1], f[2])
        t0 = fp_sub(t1, t0)
        t1 = fp_add(t2, t3)
        t2 = fp_mul(t0, g[0])
        t3 = fp_mul(t0, g[1])
        t4 = fp_mul(t1, g[0])
        t3 = fp_add(t3, t4)
        return [fp_mul(g[0], a), fp_mul(g[1], a), fp_sub(0, t2), fp_sub(0, t3)], fp_sqr(a)
    
    else:

        m = n - (n // 2)    # Ceiling of n/2
        
        # Recursively call to reciprocal
        g, a = reciprocal(f, flen, m)

        """
        # Basic idea
        t = poly_mul_modxn(n, f[:n], n, g, m)                                   #  f * g          mod x^n
        t = [fp_sub(t[0], a) ] + t[1:n]                                         # (f * g - a)     mod x^n
        assert( [0]*m == t[:m] )                                                # the first m coefficients are equal zero
        t = poly_mul_modxn(n - m, g, m, t[m:], n - m)                           # (f * g - a) * g mod x^n

        # Finally, the reciprocal is a*g - (f*g - a)*g mod x^n
        t = [ fp_mul(a, g[i]) for i in range(0, m, 1) ] + [ fp_sub(0, t[i]) for i in range(0, n - m, 1) ]
        """

        # Basic idea but saving multiplication because of f*g - a is multiple of x^m
        t = poly_mul_middle(g, m, f[:n], n)                                     #  f * g          mod x^n (the last 'm >= n - m' coefficients)
        t = poly_mul_modxn(n - m, g, m, t[(2*m - n):], n - m)                   # (f * g - a) * g mod x^n

        # Finally, the reciprocal is a*g - (f*g - a)*g mod x^n
        t = [ fp_mul(a, g[i]) for i in range(0, m, 1) ] + [ fp_sub(0, t[i]) for i in range(0, n - m, 1) ]
        
        return t, fp_sqr(a)

# Modular reduction in fp[x] with precomputation
def poly_redc(h, hlen, f):

    flen = f['deg'] + 1

    if hlen < flen:
        # Base case, h(x) mod f(x) = h(x)
        return list(h) + [0]*(flen - hlen - 1)

    elif flen == 2 and hlen == 2:

        t0 = fp_mul(h[0], f['poly'][1])     # h0 * f1
        t1 = fp_mul(h[1], f['poly'][0])     # h1 * f0
        return [ fp_sub(t0, t1) ]           # f1 * (h0 + h1*x) mod (f0 + f1*x)

    elif flen == 2 and hlen == 3:

        f0_squared = fp_sqr(f['poly'][0])       # f0^2
        f1_squared = fp_sqr(f['poly'][1])       # f1^2
        t = fp_sub(f['poly'][0], f['poly'][1])  # f0 - f1
        
        t = fp_sqr(t)                           # (f0 - f1)^2
        t = fp_sub(t, f0_squared)               # (f0 - f1)^2 - f0^2
        t = fp_sub(t, f1_squared)               # (f0 - f1)^2 - f0^2 - f1^2 = -2*f0*f1

        f0_squared = fp_add(f0_squared, f0_squared) # 2*(f0^2)
        f1_squared = fp_add(f1_squared, f1_squared) # 2*(f1^2)
        
        t0 = fp_mul(f0_squared, h[2])           # [2*(f0^2)] * h2
        t1 = fp_mul(f1_squared, h[0])           # [2*(f1^2)] * h0
        t2 = fp_mul(t, h[1])                    # [-2*f0*f1] * h1
        return [ fp_add(t0, fp_add(t1,t2))  ]   # [2 * (f1^2)] * (h0 + h1*x + h2*x^2) mod (f0 + f1*x)
    
    elif flen == 3 and hlen == 3:

        f2h1 = fp_mul(f['poly'][2], h[1])
        f2h0 = fp_mul(f['poly'][2], h[0])
        f1h2 = fp_mul(f['poly'][1], h[2])
        f0h2 = fp_mul(f['poly'][0], h[2])
        return [fp_sub(f2h0, f0h2), fp_sub(f2h1, f1h2)]
        """ 
    elif flen == hlen:

        t0 = [ fp_mul(h[inner], f['poly'][fdeg - 1]) for inner in range(0, fdeg - 1, 1) ]
        t1 = [ fp_mul(h[fdeg - 1], f['poly'][inner]) for inner in range(0, fdeg - 1, 1) ]
        return [ fp_sub(t1[inner], t0[inner]) for inner in range(0, fdeg - 1, 1) ]
        """
    else:

        H = [ h[i] for i in range(hlen - 1, -1, -1) ]                       # x^deg(h) * h(x^-1)
        Q = list(f['reciprocal'])                                           # (1/F) mod x^(deg(f) - deg(h))

        # (H/F) mod x^(deg(f) - deg(h))
        Q = poly_mul_modxn(hlen - flen + 1, Q[:(hlen - flen + 1)], hlen - flen + 1, H[:(hlen - flen + 1)], hlen - flen + 1)

        q = [ Q[i] for i in range(hlen - flen, -1 ,-1) ]                    # x^deg(Q) * Q(x^-1) is the quotient 
        qf= poly_mul_modxn(flen - 1, q, hlen - flen + 1, f['poly'], flen)   # (q * f) (x) has degree equals deg(h)
        
        # a*h(x) - (q * f)(x) will gives a polynomial of degree (deg(f) - 1)
        return [ fp_sub(fp_mul(f['a'], h[i]), qf[i]) for i in range(0, flen - 1, 1) ]

# Reciprocal tree of a given product tree
def reciprocal_tree(r, glen, ptree_f, n):

    if (n == 0):
        # Super special base case (nothing to do)
        return {'left':None, 'right':None, 'poly':[1], 'deg': 0, 'reciprocal':[1], 'a':1}

    if (n == 1):
        # We are in the leaf
        return {'left':None,\
                'right':None,\
                'poly':list(ptree_f['poly']),\
                'deg':ptree_f['deg'],\
                'reciprocal':[1],\
                'a':1}

    if ptree_f['deg'] == 2 and (glen == 3):
        # Not required because of our redcution base cases
        return {'left':ptree_f['left'],\
                'right':ptree_f['right'],\
                'poly':list(ptree_f['poly']),\
                'deg':ptree_f['deg'],\
                'reciprocal':[1],\
                'a':1}

    # At this point, we deal with a general case
    flen = ptree_f['deg'] + 1
    #R, A = reciprocal(ptree_f['poly'][::-1], flen, glen - flen + 1  )
    if r['rdeg'] <= (glen - flen):

        # This case requires a reciprocal computation
        R, A = reciprocal(ptree_f['poly'][::-1], flen, glen - flen + 1  )

    else:

        # This case allows to use a polynomial multiplication modulo x^{glen - flen + 1}
        A = r['a']
        R = poly_mul_modxn(glen - flen + 1, r['rpoly'], r['rdeg'] + 1, r['fpoly'][::-1], r['fdeg'] + 1 )

    # Now, we proceed by recusion calls
    m = n - (n//2)
    left = reciprocal_tree({'rpoly':R, 'rdeg':(glen - flen), 'fpoly':ptree_f['right']['poly'], 'fdeg':ptree_f['right']['deg'], 'a':A}, flen - 1, ptree_f['left'], m)
    right= reciprocal_tree({'rpoly':R, 'rdeg':(glen - flen), 'fpoly':ptree_f['left']['poly'], 'fdeg':ptree_f['left']['deg'], 'a':A}, flen - 1, ptree_f['right'], n - m)

    return {'left':left,\
            'right':right,\
            'poly':list(ptree_f['poly']),\
            'deg':ptree_f['deg'],\
            'reciprocal':R,\
            'a':A}

# Next function computes g(x) mod f_1(x), ..., g(x) mod f_n(x)
def multieval_unscaled(g, glen, ptree_f, n):

    if n == 0:
        return [ [1] ]

    g_mod = poly_redc(g, glen, ptree_f)

    if n == 1:

        # Now, we have g corresponds with the initial G but now it is modulus a leaf of the product tree of f
        return [g_mod]

    else:

        m = n - (n // 2)
        # Reducing g(x) modulo the current node polynomial
        left = multieval_unscaled(g_mod, ptree_f['deg'], ptree_f['left'], m )
        right = multieval_unscaled(g_mod, ptree_f['deg'], ptree_f['right'], n - m )
        return left + right 

# Next functions computes the scaled remainder tree
def multieval_scaled(g, glen, f, flen, ptree_f, n):

    if n == 0:
        return [ [1] ]

    #fg = poly_mul_middle(f, flen, g, glen)
    if flen == n and glen == n and n > 1:
        fg = list(g)
    else:
        fg = poly_mul_middle(f, flen, g, glen)

    if n == 1:
        # The last coefficient should be the desire modular reduction with linear modulus
        if fg != []:
            return [ [fg[-1]] ]
        else:
            return [ [1] ]

    m = n - (n//2)
    left = multieval_scaled(fg, flen, ptree_f['right']['poly'], ptree_f['right']['deg'] + 1, ptree_f['left'] , m)
    right= multieval_scaled(fg, flen, ptree_f['left']['poly'] ,  ptree_f['left']['deg'] + 1, ptree_f['right'], n - m)
    return left + right
