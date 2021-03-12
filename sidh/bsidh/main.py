from sidh.framework import *
from sidh.bsidh.strategy import *

if setting.verbose:
    verb = '-suitable'
else:
    verb = '-classical'

try:

    f = open(
        strategy_data
        + setting.algorithm
        + '-'
        + setting.prime
        + '-'
        + setting.formulaes
        + verb
    )
    print("// Strategies to be read from a file")

    # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to include case l=2 and l=4]
    tmp = f.readline()
    tmp = [int(b) for b in tmp.split()]
    Sp = list(tmp)
    # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    tmp = f.readline()
    tmp = [int(b) for b in tmp.split()]
    Sm = list(tmp)

    f.close()

except IOError:

    print("// Strategies to be computed")
    # List of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to include case l=2 and l=4]
    Sp, Cp = dynamic_programming_algorithm(SIDp[::-1], len(SIDp))

    # List of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    Sm, Cm = dynamic_programming_algorithm(SIDm[::-1], len(SIDm))

    f = open(
        strategy_data
        + setting.algorithm
        + '-'
        + setting.prime
        + '-'
        + setting.formulaes
        + verb,
        'w',
    )

    f.writelines(' '.join([str(tmp) for tmp in Sp]) + '\n')
    f.writelines(' '.join([str(tmp) for tmp in Sm]) + '\n')

    f.close()

print(
    "// All the experiments are assuming S = %1.6f x M and a = %1.6f x M. The measures are given in millions of field operations.\n"
    % (SQR, ADD)
)


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

print("public_coeff := 0x%X + i * 0x%X;\n" % (a[0], a[1]))

print(
    "// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
)
print("// Public Key Generation")
print(
    "// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
)

# Alice's side
print("// Private key corresponding to Alice")
a_private = random_key(p - 1)
set_zero_ops()
Ra = Ladder3pt(a_private, PA, QA, PQA, A)
print("// sk_a := 0x%X;" % a_private)
RUNNING_TIME = get_ops()
print(
    "// kernel point generator cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;\n"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)
print("// Public key corresponding to Alice")
set_zero_ops()
a_public, PB_a, QB_a, PQB_a = evaluate_strategy(
    True, PB, QB, PQB, A, Ra, SIDp[::-1], Sp, len(SIDp)
)
RUNNING_TIME = get_ops()
print(
    "// isogeny evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)

a_curve = coeff(a_public)
print("pk_a := 0x%X + i * 0x%X;\n" % (a_curve[0], a_curve[1]))
# print("B := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (a_curve[0], a_curve[1]) )
# print("assert(Random(B) * (p + 1) eq B!0);")

print("// Private key corresponding to Bob")
b_private = random_key(p - 1)
print("// sk_b := 0x%X;" % b_private)
set_zero_ops()
Rb = Ladder3pt(b_private, PB, QB, PQB, A)
RUNNING_TIME = get_ops()
print(
    "// kernel point generator cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;\n"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)
print("// Public key corresponding to Bob")
set_zero_ops()
b_public, PA_b, QA_b, PQA_b = evaluate_strategy(
    True, PA, QA, PQA, A, Rb, SIDm[::-1], Sm, len(SIDm)
)
RUNNING_TIME = get_ops()
print(
    "// isogeny evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)

b_curve = coeff(b_public)
print("pk_b := 0x%X + i * 0x%X;\n" % (b_curve[0], b_curve[1]))
# print("B := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (b_curve[0], b_curve[1]) )
# print("assert(Random(B) * (p + 1) eq B!0);")

# ======================================================
print(
    "\n// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
)
print("// Secret Sharing Computation")
print(
    "// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
)

print("// Secret sharing corresponding to Alice")
set_zero_ops()
RB_a = Ladder3pt(a_private, PA_b, QA_b, PQA_b, b_public)
RUNNING_TIME = get_ops()
print(
    "// kernel point generator cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)
set_zero_ops()
ss_a, _, _, _ = evaluate_strategy(
    False, PB, QB, PQB, b_public, RB_a, SIDp[::-1], Sp, len(SIDp)
)
RUNNING_TIME = get_ops()
print(
    "// isogeny evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)

ss_a_curve = coeff(ss_a)
print("ss_a := 0x%X + i * 0x%X;\n" % (ss_a_curve[0], ss_a_curve[1]))
# ss_ja = jinvariant(ss_a)
# print("B := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (ss_a_curve[0], ss_a_curve[1]) )
# print("jB := 0x%X + i * 0x%X;" % (ss_ja[0], ss_ja[1]) )
# print("assert(Random(B) * (p + 1) eq B!0);")
# print("assert(jInvariant(B) eq jB);")

print("// Secret sharing corresponding to Bob")
set_zero_ops()
RA_b = Ladder3pt(b_private, PB_a, QB_a, PQB_a, a_public)
RUNNING_TIME = get_ops()
print(
    "// kernel point generator cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)
set_zero_ops()
ss_b, _, _, _ = evaluate_strategy(
    False, PA, QA, PQA, a_public, RA_b, SIDm[::-1], Sm, len(SIDm)
)
RUNNING_TIME = get_ops()
print(
    "// isogeny evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
    % (
        RUNNING_TIME[0] / (10.0 ** 6),
        RUNNING_TIME[1] / (10.0 ** 6),
        RUNNING_TIME[2] / (10.0 ** 6),
        measure(RUNNING_TIME) / (10.0 ** 6),
    )
)

ss_b_curve = coeff(ss_b)
print("ss_b := 0x%X + i * 0x%X;\n" % (ss_b_curve[0], ss_b_curve[1]))
# ss_jb = jinvariant(ss_b)
# print("B := EllipticCurve(x^3 + (0x%X + i * 0x%X )* x^2 + x);" % (ss_b_curve[0], ss_b_curve[1]) )
# print("jB := 0x%X + i * 0x%X;" % (ss_jb[0], ss_jb[1]) )
# print("assert(Random(B) * (p + 1) eq B!0);")
# print("assert(jInvariant(B) eq jB);")

print(
    "\n// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
)
if ss_a_curve == ss_b_curve:
    print('\x1b[0;30;43m' + '\"Successfully passed!\";' + '\x1b[0m')
else:
    print(
        '\x1b[0;30;41m'
        + '\"Great Scott!... The sky is falling. NOT PASSED!!!\"'
        + '\x1b[0m'
    )
