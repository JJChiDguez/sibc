from sidh.framework import *

if setting.formulaes == 'tvelu':
    from sidh.bsidh.tvelu import *

elif setting.formulaes == 'svelu':
    from sidh.bsidh.svelu import *

else:
    from sidh.bsidh.hvelu import *

print("p := 0x%X;" % p)
print("fp := GF(p);")
print("_<x> := PolynomialRing(fp);")
print("fp2<i> := ext<fp | x^2 + 1>;")
print("Pr<x> := PolynomialRing(fp2);")

# Reading public generators points
f = open(gen_data + setting.prime)

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

print("E := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (a[0], a[1]))

S = [list(PA[0]), list(PA[1])]
T = [list(QA[0]), list(QA[1])]
ST = [list(PQA[0]), list(PQA[1])]

for i in range(0, np, 1):
    for idx in range(0, Ep[i], 1):
        S = xMUL(S, A, i)
        T = xMUL(T, A, i)
        ST = xMUL(ST, A, i)

assert isinfinity(S)
assert isinfinity(T)
assert isinfinity(ST)

print("\n// Verifying torsion-(p + 1) points")
print("// x([p + 1]PA) = (1:0)?\t", isinfinity(S))
print("// x([p + 1]QA) = (1:0)?\t", isinfinity(T))
print("// x([p + 1]PQA) = (1:0)?\t", isinfinity(ST))

S = [list(PB[0]), list(PB[1])]
T = [list(QB[0]), list(QB[1])]
ST = [list(PQB[0]), list(PQB[1])]

for i in range(np, np + nm, 1):
    for idx in range(0, Em[i - np], 1):
        S = xMUL(S, A, i)
        T = xMUL(T, A, i)
        ST = xMUL(ST, A, i)

assert isinfinity(S)
assert isinfinity(T)
assert isinfinity(ST)
print("\n// Verifying torsion-(p - 1) points")
print("// x([p - 1]PB) = (1:0)?\t", isinfinity(S))
print("// x([p - 1]QB) = (1:0)?\t", isinfinity(T))
print("// x([p - 1]PQB) = (1:0)?\t", isinfinity(ST))

# Case (p + 1)
S = [list(PA[0]), list(PA[1])]
T = [list(QA[0]), list(QA[1])]
ST = [list(PQA[0]), list(PQA[1])]

assert isinfinity(S) == False
assert isinfinity(T) == False
assert isinfinity(ST) == False

for i in range(0, np, 1):
    for idx in range(0, Ep[i] - 1, 1):
        S = xMUL(S, A, i)
        T = xMUL(T, A, i)
        ST = xMUL(ST, A, i)

print("\n// Verifying orders")
assert isfull_order(prime_factors(S, A, range(0, np, 1)))
assert isfull_order(prime_factors(T, A, range(0, np, 1)))
assert isfull_order(prime_factors(ST, A, range(0, np, 1)))
print(
    "// PA is a full order point?\t",
    isfull_order(prime_factors(S, A, range(0, np, 1))),
)
print(
    "// QA is a full order point?\t",
    isfull_order(prime_factors(T, A, range(0, np, 1))),
)
print(
    "// QPA is a full order point?\t",
    isfull_order(prime_factors(ST, A, range(0, np, 1))),
)

# Case (p - 1)
S = [list(PB[0]), list(PB[1])]
T = [list(QB[0]), list(QB[1])]
ST = [list(PQB[0]), list(PQB[1])]

assert isinfinity(S) == False
assert isinfinity(T) == False
assert isinfinity(ST) == False

for i in range(np, np + nm, 1):
    for idx in range(0, Em[i - np] - 1, 1):
        S = xMUL(S, A, i)
        T = xMUL(T, A, i)
        ST = xMUL(ST, A, i)

print("\n// Verifying orders")
assert isfull_order(prime_factors(S, A, range(np, np + nm, 1)))
assert isfull_order(prime_factors(T, A, range(np, np + nm, 1)))
assert isfull_order(prime_factors(ST, A, range(np, np + nm, 1)))
print(
    "// PB is a full order point?\t",
    isfull_order(prime_factors(S, A, range(np, np + nm, 1))),
)
print(
    "// QB is a full order point?\t",
    isfull_order(prime_factors(T, A, range(np, np + nm, 1))),
)
print(
    "// QPB is a full order point?\t",
    isfull_order(prime_factors(ST, A, range(np, np + nm, 1))),
)

# --------------------------------------------------------------------------------------------------------------------
# Three point ladder: case (p + 1)
S = [list(PA[0]), list(PA[1])]
T = [list(QA[0]), list(QA[1])]
ST = [list(PQA[0]), list(PQA[1])]

assert isinfinity(S) == False
assert isinfinity(T) == False

for i in range(0, np, 1):
    for idx in range(0, Ep[i] - 1, 1):
        S = xMUL(S, A, i)
        T = xMUL(T, A, i)
        ST = xMUL(ST, A, i)

k = random.randint(0, p)
R = Ladder3pt(k, S, T, ST, A)
# print("k := 0x%X;" % k)
# print("boolR, R := IsPoint(E, (0x%X + i * 0x%X) / (0x%X + i * 0x%X));" % (R[0][0], R[0][1], R[1][0], R[1][1]))

T_p = [list(R[0]), list(R[1])]
T_m = [list(S[0]), list(S[1])]
print(
    "\n// Now, we proceed by performing xISOG with input curve equals the output curve of the previous one experiment, and using torsion-(p + 1) points"
)
for idx in range(0, np, 1):

    # -------------------------------------------------------------
    # Random kernel point
    Tp = list(T_p)
    for i in range(idx + 1, np, 1):
        Tp = xMUL(Tp, A, i)

    print("// l:\t%7d |" % global_L[idx], end="")
    total_cost = [0, 0, 0]

    if setting.formulaes != 'tvelu':

        if setting.verbose:
            set_parameters_velu(sJ_list[idx], sI_list[idx], idx)

        else:
            # -------------------------------------------------------------
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in KPs, xISOG, and xEVAL
            if global_L[idx] <= 4:
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
    Tp[0], A[0] = fp2_cswap(Tp[0], A[0], global_L[idx] == 4)
    Tp[1], A[1] = fp2_cswap(Tp[1], A[1], global_L[idx] == 4)
    B = xISOG(A, idx)
    Tp[0], A[0] = fp2_cswap(Tp[0], A[0], global_L[idx] == 4)
    Tp[1], A[1] = fp2_cswap(Tp[1], A[1], global_L[idx] == 4)
    show_ops("xISOG", 1.0, 0.0, False)
    t = get_ops()
    total_cost[0] += t[0]
    total_cost[1] += t[1]
    total_cost[2] += t[2]

    # -------------------------------------------------------------
    # xEVAL: kernel point determined by the next isogeny evaluation
    set_zero_ops()
    if (
        setting.formulaes == 'tvelu'
        or (setting.formulaes == 'hvelu' and global_L[idx] <= HYBRID_BOUND)
        or (global_L[idx] == 4)
    ):
        T_p = xEVAL(T_p, idx)
    else:
        T_p = xEVAL(T_p, A)

    # xEVAL bench
    set_zero_ops()
    if (
        setting.formulaes == 'tvelu'
        or (setting.formulaes == 'hvelu' and global_L[idx] <= HYBRID_BOUND)
        or (global_L[idx] == 4)
    ):
        T_m = xEVAL(T_m, idx)
    else:
        T_m = xEVAL(T_m, A)

    show_ops("xEVAL", 1.0, 0.0, False)
    t = get_ops()
    total_cost[0] += t[0]
    total_cost[1] += t[1]
    total_cost[2] += t[2]
    print("|| cost: %8d" % (total_cost[0] + total_cost[1]), end=" ")
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
# print("B := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(A))
a = coeff(A)
print("B := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (a[0], a[1]))
print("assert(Random(B) * (p + 1) eq B!0);")

# --------------------------------------------------------------------------------------------------------------------
A = [[0x8, 0x0], [0x4, 0x0]]
# Three point ladder: case (p - 1)
S = [list(PB[0]), list(PB[1])]
T = [list(QB[0]), list(QB[1])]
ST = [list(PQB[0]), list(PQB[1])]

assert isinfinity(S) == False
assert isinfinity(T) == False

for i in range(np, np + nm, 1):
    for idx in range(0, Em[i - np] - 1, 1):
        S = xMUL(S, A, i)
        T = xMUL(T, A, i)
        ST = xMUL(ST, A, i)

k = random.randint(0, p)
R = Ladder3pt(k, S, T, ST, A)
# print("k := 0x%X;" % k)
# print("boolR, R := IsPoint(E, (0x%X + i * 0x%X) / (0x%X + i * 0x%X));" % (R[0][0], R[0][1], R[1][0], R[1][1]))

T_p = [list(R[0]), list(R[1])]
T_m = [list(S[0]), list(S[1])]
print(
    "\n// Now, we proceed by performing xISOG with input curve equals the output curve of the previous one experiment, and using torsion-(p - 1) points"
)
for idx in range(np, np + nm, 1):

    # -------------------------------------------------------------
    # Random kernel point
    Tp = list(T_p)
    for i in range(idx + 1, np + nm, 1):
        Tp = xMUL(Tp, A, i)

    print("// l:\t%7d |" % global_L[idx], end="")
    total_cost = [0, 0, 0]

    if setting.formulaes != 'tvelu':

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
    if setting.formulaes == 'tvelu' or (
        setting.formulaes == 'hvelu' and global_L[idx] <= HYBRID_BOUND
    ):
        T_p = xEVAL(T_p, idx)
    else:
        T_p = xEVAL(T_p, A)

    # xEVAL bench
    set_zero_ops()
    if setting.formulaes == 'tvelu' or (
        setting.formulaes == 'hvelu' and global_L[idx] <= HYBRID_BOUND
    ):
        T_m = xEVAL(T_m, idx)
    else:
        T_m = xEVAL(T_m, A)

    show_ops("xEVAL", 1.0, 0.0, False)
    t = get_ops()
    total_cost[0] += t[0]
    total_cost[1] += t[1]
    total_cost[2] += t[2]
    print("|| cost: %8d" % (total_cost[0] + total_cost[1]), end=" ")
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
a = coeff(A)
print("B := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (a[0], a[1]))
print("assert(Random(B) * (p + 1) eq B!0);")
# """

print(
    "\n\"If no errors were showed using magma calculator, then all experiments were successful passed!\";"
)
print("// copy and paste it at http://magma.maths.usyd.edu.au/calc/\n")

"""
from bsidh.poly_redc import *
flen = random.randint(1, 32)
glen = random.randint(1, 32)
f = [ [random.randint(0, p - 1), random.randint(0, p - 1)]  for i in range(0, flen, 1) ]
g = [ [random.randint(0, p - 1), random.randint(0, p - 1)]  for i in range(0, glen, 1) ]
fg = poly_mul(f, flen, g, glen)
print("f := ", f, ";")
print("g := ", g, ";")
print("fg := ", fg, ";")
print("Pr![ fp2!fi : fi in f] * Pr![ fp2!gi : gi in g] eq Pr![ fp2!fgi : fgi in fg];");
"""
