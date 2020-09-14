from bsidh.strategy import *

if setting.verbose:
    verb = '-suitable'
else:
    verb = '-classical'

try:

    f = open('./strategies/' + setting.algorithm + '-' + setting.prime  + '-' + setting.formulaes + verb)
    print("// Strategies to be read from a file")

    # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to include case l=2 and l=4]
    tmp = f.readline()
    tmp = [ int(b) for b in tmp.split() ]
    Sp = list(tmp)
    # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    tmp = f.readline()
    tmp = [ int(b) for b in tmp.split() ]
    Sm = list(tmp)

    f.close()

except IOError:

    print("// Strategies to be computed")
    # List of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to include case l=2 and l=4]
    Sp, Cp = dynamic_programming_algorithm( SIDp[::-1], len(SIDp))
    
    # List of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    Sm, Cm = dynamic_programming_algorithm( SIDm[::-1], len(SIDm))

    f = open('./strategies/' + setting.algorithm + '-' + setting.prime  + '-' + setting.formulaes + verb,'w')

    f.writelines(' '.join([ str(tmp) for tmp in Sp]) + '\n')
    f.writelines(' '.join([ str(tmp) for tmp in Sm]) + '\n')

    f.close()

# List of strategies
S_out = [list(Sp), list(Sm)]
filename = './figures/' + setting.algorithm + '-' + setting.prime  + '-' + setting.formulaes + verb
