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

# ---
k = 3   # Number of rows (format of the list)
# ---

print("#ifndef _IJK_%s_H_" % setting.prime)
print("#define _IJK_%s_H_" % setting.prime)
print("")
assert(n == len(L))

print("#ifdef _MONT_C_CODE_")
print("// The list of the bitlength of each SOP")
printl("static uint64_t bL[]", [bitlength(l) for l in L], n // k + 1)
print("#endif")

print("")
print("#ifdef _ISOG_H_")
print("\n// The list of Small Odd Primes (SOPs) is stored such that l_0 < l_1 < ... < l_{n-1}")
printl("static uint64_t L[]", L, n // k + 1)

assert(n == len(sI_list))
assert(n == len(sJ_list))

sK_list = []
for i in range(0, n, 1):
    assert(sJ_list[i] <= sI_list[i])
    sK_list.append(((L[i] - 2 - 4*sJ_list[i]*sI_list[i] - 1) // 2) + 1)
    assert(sK_list[i] >= 0)

print("")
print("#ifndef _C_CODE_")
print("// Sizes for the sets I, J, and K required in the new velusqrt formulae")
printl("static int sizeI[]", sI_list, n // k + 1)
printl("static int sizeJ[]", sJ_list, n // k + 1)
printl("static int sizeK[]", sK_list, n // k + 1)
print("#endif")

print("")
print("#define sI_max %d" % (max(sI_list)))
print("#define sJ_max %d" % (max(sJ_list)))
print("#define sK_max %d" % (L[n-1] // 2 + 1))
print("#endif")

print("\n#endif /* required framework for the #I, #J, and #K, which is used in new velusqrt fomurlae on CSIDH-%s */" % setting.prime[1:] )
