from sys import argv
from sympy import Symbol, symbols, sympify, expand
from sympy import Mul, Add, Pow
from itertools import combinations

symb_dict = {} # entries are (varname, symbolname), presently assumed equal
symb_dict['Freq'] = ('w', 'w')
symb_dict['Ener'] = ('E', 'E') # fixme: has to be single character
symb_dict['Fock'] = ('F', 'F')
symb_dict['Ovlp'] = ('S', 'S')
symb_dict['Dens'] = ('D', 'D')
symb_dict['Oidt'] = ('iSt', 'iSt') # Notation: iSt represents i*dS/dt
symb_dict['Didt'] = ('iDt', 'iDt') #

expr_dict = {} # entries are (variable name, sympy expression)
expr_dict['Pulay'] = ('Wt', 'D*F*D + 0.5*iDt*S*D - 0.5*D*S*iDt')
expr_dict['Scfeq'] = ('Yt', 'F*D*S - S*D*F - ' + 
                           'S*iDt*S - 0.5*iSt*D*S - 0.5*S*D*iSt')
expr_dict['Idemp'] = ('Z', 'D*S*D - D')
expr_dict['ScfeqMul'] = ('ya', 'Da*S*D - D*S*Da')
expr_dict['IdempMul'] = ('zat', 'Fa*D*S + S*D*Fa - Fa - F*D*Sa - Sa*D*F + ' +
                         'S*iDt*Sa - Sa*iDt*S + 0.5*iSt*D*Sa - 0.5*Sa*D*iSt')

def var_name(key):
    '''Return variable name for given key in symb_dict or expr_dict'''
    if key in symb_dict: return symb_dict[key][0]
    if key in expr_dict: return expr_dict[key][0]
    return ''

lagr_dict = {} # terms in the Lagrangian
lagr_dict['Energy'] = var_name('Ener') + '0at'
lagr_dict['Pulay'] = '%sa*%s' % (var_name('Ovlp'), var_name('Pulay'))
lagr_dict['Idempotency'] = '%s*%s' % (var_name('IdempMul'), var_name('Idemp'))
lagr_dict['Scf Eqn'] = '%s*%s' % (var_name('ScfeqMul'), var_name('Scfeq'))

def pert_name(name, pseq):
    '''Return list of perturbed names for perturbations in pseq'''
    return [name + ''.join(map(str,p)) for p in pseq]

def comb_k(seq, k):
    """Return all combinations of seq up to order k,
    where the first element of seq is always present"""
    kcomb = []
    for order in range(k):
        for comb in combinations(seq[1:], order):
            kcomb.append((seq[0],)+comb)
    return kcomb

def comb_kn(seq, k, n):
    """Return all combinations of seq up to order k in the first
    element and up to order n in the following elements.  If seq is
    ordered so is the return vale."""
    kncomb = []
    for order in range(n+1):
        ncomb = []
        for comb in combinations(seq[1:], order):
            ncomb.append(comb)
        kncomb.extend(ncomb)
        if order < k:
            for comb in ncomb:
                kncomb.append((seq[0],)+comb)
    return kncomb

def assign_symvars(vseq, pseq, is_commutative=True):
    '''Return string to use with "exec" for dynamic assignment of
    variables to sympy symbols.  "exec exe_str" will assign the sympy
    symbols in pseq to the python variables in vseq.'''
    exe_str = ''
    for vname, sname in zip(vseq,pseq):
        cstr = str(is_commutative)
        estr = "%s = Symbol('%s',commutative=%s)" % (vname, sname, cstr)
        exe_str = exe_str + estr + "; "
    return exe_str

# Set total order m and type of n+k+1 rule
m, k, n = map(int,argv[1:4])
print 'Input: (m, k, n) =', (m, k, n), '\n' # fixme: do proper input processing
assert(m == k+n+1), "k,n out of range (m, k, n) = (%d, %d, %d)" % (m, k, n)
assert(k <= m/2  ), "k,n out of range (m, k, n) = (%d, %d, %d)" % (m, k, n)
assert(n <= m-1  ), "k,n out of range (m, k, n) = (%d, %d, %d)" % (m, k, n)

# Represent perturbations as symbols a, b, ...
pert = tuple([ Symbol(chr(ord('a')+i)) for i in range(m) ])

print 'Quasi-energy derivative of order', m 
print 'with respect to perturbations', pert
print 'using k,n-rule with (k, n) = (%d, %d)\n' % (k, n)

print 'Perturbed densities needed'
print 'up to order (k = %d) in %s : %s' % (
    k, pert[:1], pert_name(var_name('Dens'), comb_k(pert, k))) 
if m > 1:
    print 'up to order (n = %d) in %s : %s' % (
        n, pert[1:], pert_name(var_name('Dens'), comb_kn(pert[1:], n, n))[1:])
print

# Symbols for frequencies 
var_list = pert_name(var_name('Freq'), comb_kn(pert, k, n))
exec assign_symvars(var_list, var_list)

# Symbols for energy derivatives eg E1a = d/da(dE/dD) 
for i in range(m+1):
    ename  = var_name('Ener') + str(i)
    var_list = pert_name(ename, comb_kn(pert, m, m))
    exec assign_symvars(var_list, var_list)

# Symbols for overlap derivatives
var_list = pert_name(var_name('Ovlp'), comb_kn(pert, m, m))
exec assign_symvars(var_list, var_list, is_commutative=False)

# Symbols for derivatives appearing in Multiplier terms
for key in ['Fock', 'Dens', 'Oidt', 'Didt']:
    var_list = pert_name(var_name(key), comb_kn(pert, k, n))
    exec assign_symvars(var_list, var_list, is_commutative=False)

dt2w = {} # dict to use with sympy.subs
wlist = pert_name(var_name('Freq'), comb_kn(pert, k, n))
tlist = pert_name(var_name('Didt'), comb_kn(pert, k, n))
slist = pert_name(var_name('Dens'), comb_kn(pert, k, n))
for tstr, wstr, sstr in zip(tlist, wlist, slist):
    dt2w[eval(tstr)] = eval('%s*%s' % (wstr, sstr)) # iDta = wa*Da etc
tlist = pert_name(var_name('Oidt'), comb_kn(pert, k, n))
slist = pert_name(var_name('Ovlp'), comb_kn(pert, k, n))
for tstr, wstr, sstr in zip(tlist, wlist, slist):
    dt2w[eval(tstr)] = eval('%s*%s' % (wstr, sstr)) # iSta = wa*Sa etc
dt2w[eval(var_name('Oidt'))] = 0 # iSt=0
dt2w[eval(var_name('Didt'))] = 0 # iDt=0
# overlap derivatives that do not contribute to the multiplier terms
for knpert in comb_kn(pert, k, n):
    spert = tuple(sorted(set(pert) - set(knpert)))
    sstr  = var_name('Ovlp') + ''.join(map(str, spert))
    dt2w[eval(sstr)] = 0

def sdiff(symb, var):
    '''Differentiation of symb with respect to variable var. (both are
    assumed to be sympy symbols)

    The derivative is represented by Symbol(symb.name+var.name). If
    the new derivative symbol does not exist it is assumed to be zero:
    sdiff(Da, b) --> Dab if Symbol('Dab') exists, 0 other whise

    A digit as second character in symb.name indicates dependency on
    the density and the chain rule is applied:
    sdiff(E1a, b) --> E1ab + E2a * sdiff(D, b)'''
    e = sympify(0) # if none of the derivative symbols exist return 0
    try: # partial derivative
        e += eval(symb.name + var.name)
    except: # fixme: be more specific (this catches all errors)
        pass # eval() fails with NameError if derivate symbol is undefined
    try: # inner derivative for chain rule
        istr = symb.name[0] + str(int(symb.name[1])+1) + symb.name[2:]
        e += eval(istr)*sdiff(eval(var_name('Dens')),var)
    except: # fixme: be more specific (this catches all errors)
        pass # int() failes with ValueError if second char is not integer
    return e

def udiff(expr, var):
    '''Univariate differentiation of sympy expression expr with
    respect to sympy Symbol var'''
    if expr.is_Symbol:
        return sdiff(expr, var)
    elif expr.is_Number:
        return sympify(0)
    elif expr.is_Add:
        return udiff(expr.args[0], var) + udiff(Add(*expr.args[1:]), var)
    elif expr.is_Mul:
        return ( udiff(expr.args[0], var) * Mul(*expr.args[1:]) + 
                 expr.args[0] * udiff(Mul(*expr.args[1:]), var) )
    elif expr.is_Pow:
        return ( udiff(expr.args[0], var) * 
                 Pow(expr.args[0], expr.args[1] - 1) * expr.args[1] )
    else:
        print 'udiff error: unknown expression type'

def mdiff(expr, mvar):
    '''Multivariate differentiation of expr (sympy expression) with
    respect to variables in mvar (sequence of sympy Symbols)'''
    ed = expr
    for var in mvar:
        ed = udiff(ed, var)
    return expand(ed)

print 'Lagrangian:\nLat = %s\n' % ' - '.join([lagr_dict[term]
  for term in ['Energy', 'Pulay', 'Scf Eqn', 'Idempotency']])

# Set variables in Lagrangian
nlist = ['Pulay', 'Scfeq', 'Idemp'] # Scfeq and Idemp needed for RHS set up
if k>0:
    nlist.extend(['ScfeqMul', 'IdempMul'])
for name in nlist:
    print 'Assign %s expression to variable:' % name
    var, expr = expr_dict[name]
    exec("%s=%s" % (var, expr))
    print '%s = %s' % (var, expr)
    print

print 'Set up right hand side and solve response equations'
print 'Construct perturbed %s and %s\n' % (var_name('Dens'), var_name('Fock'))
for npert in comb_kn(pert, k, n)[1:]: # [1:] skips first unperturbed 
    pstr  = ''.join(map(str,npert))
    dstr = var_name('Dens') + pstr
    fstr = var_name('Fock') + pstr
    print 'Set up response equation for: ', dstr
    print 'Initialize particular solution Dp:'
    print 'Dp =', mdiff(eval(var_name('Idemp')), npert).subs(eval(dstr),0)
    print 'Project: ', 'Dp = Dp - D*S*Dp - Dp*S*D'
    print 'Construct %s with %s = Dp' % (fstr, dstr)
    print 'Construct Right Hand Side with %s = Dp' % dstr
    print 'RHS = ', mdiff(eval(var_name('Scfeq')), npert).subs(dt2w)
    print 'Solve for homogenous solution Dh'
    print 'Add Dh contribution to complete: ', dstr, fstr
    print

# variable name to store string representation of derivative Lagrangian
# Notation: Lakbcd...n = L_{k,n}^{abcd...}
lstr = 'L%s%d' % (pert[0], k)
if m > 1:
    lstr += ''.join(map(str, pert[1:])) + str(n)

print 'Energy terms:\n'
estr = 'E%s%d' % (pert[0], k) 
if m > 1:
    estr += ''.join(map(str ,pert[1:])) + str(n)
exec "%s = '%s'" % (lstr, estr)
exec '%s = mdiff(%s, pert)' % (estr, var_name('Ener') + '0')
print '%s = %s\n' % (estr, eval(estr))

term = 'Pulay'
print '%s energy terms:\n' % term
sname, wname = lagr_dict[term].split('*')
for wpert in comb_kn(pert, k, n):
    spert = tuple(sorted(set(pert) - set(wpert))) # complement to wpert
    wstr  = wname.strip('t') + ''.join(map(str, wpert))
    sstr  = sname.strip('t') + ''.join(map(str, spert[1:]))
    exec "%s += ' - %s*%s'" % (lstr, sstr, wstr) 
    print '%s * %s' % (sstr, wstr)
    exec '%s = %s' % (wstr, mdiff(eval(wname), wpert).subs(dt2w))
    print '%s = %s\n' % (wstr, eval(wstr)) 

for term in ['Pulay', 'Scf Eqn', 'Idempotency']:
    print '%s multiplier terms:\n' % term
    mname, ename = lagr_dict[term].split('*')
    for mpert in comb_k(pert, k):
        epert = tuple(sorted(set(pert) - set(mpert)))
        mstr  = mname.strip('t') + ''.join(map(str, mpert[1:]))
        estr  = ename.strip('t') + str(n ) + ''.join(map(str, epert))
        exec "%s += ' - %s*%s'" % (lstr, mstr, estr) 
        print '%s * %s' % (mstr, estr)
        exec '%s = %s' % (mstr, mdiff(eval(mname), mpert[1:]).subs(dt2w))
        print '%s = %s' % (mstr, eval(mstr)) 
        exec '%s = %s' % (estr, mdiff(eval(ename), epert).subs(dt2w))
        print '%s = %s\n' % (estr, eval(estr))
    
print 'Lagrangian:\n%s = %s' % (lstr, eval(lstr))
