from sage.all import ZZ # cannot seem to import this from the right place
from sage.functions.other import binomial
from sage.rings.padics.factory import Qp, Zp, ZpCA
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.arith.functions import lcm
from sage.matrix.special import zero_matrix, identity_matrix, matrix
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.structure.factorization import Factorization
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from pycontrolledreduction import controlledreduction
import multiprocessing
ncpus =  multiprocessing.cpu_count()
import warnings
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob
from sage.misc.cachefunc import cached_function


def from_H1_to_H2(cp1, F1, tensor=False):
    """
    Input:
        cp1 - characteristic polynomial of Frob acting on H^1
        F1 - Frob approximation of the action H^1
        tensor - replace H^2 by H^1 otimes H^1
    Output:
        cp2 - characteristic polynomial of Frob acting on H^2
        F2 - Frob approximation of the action H^2

    """
    if tensor:
        return from_H1_to_H1otimesH1(cp1, F1)
    # deal with Frob approximation
    wedge_basis = [(i,j) for i in range(F1.nrows()) for j in range(i+1, F1.nrows())]
    M = FiniteRankFreeModule(F1.base_ring(), F1.nrows())
    M.basis('e')
    F1e = [M(elt) for elt in F1.columns()]
    F2e = [F1e[i].wedge(F1e[j]) for i, j in wedge_basis]
    F2 = matrix([[elt.comp()[i,j] for elt in F2e] for i,j in wedge_basis])
    # deduce the correct characteristic polynomial
    assert cp1[0] == 1, "The constant term of the characteristic polynomial must be 1"
    cp2 = Lpoly_H2(cp1)
    return cp2, F2


def tensor_charpoly(f, g):
    r"""
    INPUT:

    - ``f`` -- the characteristic polynomial of a linear transformation

    - ``g`` -- the characteristic polynomial of a linear transformation

    OUTPUT: the characteristic polynomial of the tensor product of the linear transformations

    EXAMPLES::

    sage: from crystalline_obstruction.main import tensor_charpoly
    sage: x = PolynomialRing(ZZ,"x").gen();
    sage: tensor_charpoly((x - 3) * (x + 2),  (x - 7) * (x + 5))
    x^4 - 2*x^3 - 479*x^2 - 420*x + 44100
    sage: (x - 21) * (x - 10) * (x + 14) * (x + 15)
    x^4 - 2*x^3 - 479*x^2 - 420*x + 44100

    """

    R = PolynomialRing(g.parent(), "y");
    y = R.gen();
    #x = g.parent().gen()
    A = f(y)
    B = R(g.homogenize(y))
    return B.resultant(A)


def base_change(Lpoly, r):
    R = Lpoly.parent()
    T = R.gen()
    S = PolynomialRing(R, 'u')
    u = S.gen()
    return R(Lpoly(u).resultant(u**r - T))

def Lpoly_H1otimesH1(Lpoly):
    return tensor_charpoly(Lpoly, Lpoly)

def Lpoly_H2(Lpoly):
    return ((tensor_charpoly(Lpoly, Lpoly)//base_change(Lpoly, 2))).sqrt()

def from_H1_to_H1otimesH1(cp1, F1):
    """
    Input:
        cp1 - characteristic polynomial of Frob acting on H^1
        F1 - Frob approximation of the action H^1
    Output:
        cp2 - characteristic polynomial of Frob acting on H^2
        F2 - Frob approximation of the action H^2

    """
    # deal with Frob approximation
    basis = [(i,j) for i in range(F1.nrows()) for j in range(F1.nrows())]
    genus = F1.nrows()//2
    weight = lambda pair: sum(-1 if elt < genus else 0 for elt in pair)
    basis.sort(key=weight)
    M = FiniteRankFreeModule(F1.base_ring(), F1.nrows())
    M.basis('e')
    F1e = [M(elt) for elt in F1.columns()]
    F2e = [F1e[i]*F1e[j] for i, j in basis]
    F2 = matrix([[elt.comp()[i,j] for elt in F2e] for i,j in basis])
    # deduce the correct characteristic polynomial
    assert cp1[0] == 1, "The constant term of the characteristic polynomial must be 1"
    cp2 = Lpoly_H1otimesH1(cp1)
    return cp2, F2

def find_monic_and_odd_model(f, p):
    """
    Given f in Qp[x],  return g in Qp[x] such that
    y^2 = f and y^2 = g are isomorphic over Qpbar
    and g is of odd degree and monic
    """
    if f.degree() % 2 == 0:
        # try to find weirstrass point
        roots = f.roots()
        if len(roots) == 0:
            raise NotImplementedError("Hyperelliptic curve doesn't have a weirstrass point over Qp")
        root = roots[0][0]
        x = f.variables()[0]
        f = f(x + root).reverse()
        assert f.degree() % 2 == 1
    return f.monic()


def compute_frob_matrix_and_cp_H2(f, p, prec, **kwargs):
    """
    Return a p-adic matrix approximating the action of Frob on H^2 of a surface or abelian variety,
    and its characteristic polynomial over ZZ
    Input:
        - f defining the curve or surface
        - p, prime
        - prec, a lower bound for the desired precision to run the computations, this increases the time exponentially
        - kwargs, keyword arguments to be passed to controlledreduction

    Output:
        - `prec`, the minimum digits absolute precision for approximation of the Frobenius
        - a matrix representing an approximation of Frob matrix with at least `prec` digits of absolute precision
        - characteristic polynomial of Frob on H^2
        - the number of classes omitted by working with primitive cohomology

    Note: if given two or one univariate polynomial, we will try to change the model over Qpbar,
    in order to work with an odd and monic model
    """
    K = Qp(p, prec=prec+10)
    Rf = f.parent()
    R = f.base_ring()
    if len(Rf.gens()) == 2:
        if min(f.degrees()) != 2:
            raise NotImplementedError("Affine curves must be hyperelliptic")
        x, y = f.variables()
        if f.degree(x) == 2:
            f = f.substitute(x=y, y=x)
        # Get Weierstrass equation
        # y^2 + a*y  + b == 0
        b, a, _ = map(R['x'], R['x']['y'](f).monic())
        # y^2 + a*y  + b == 0 --> (2y + a)^2 = a^2 - 4 b
        f = a**2 - 4*b
        f = find_monic_and_odd_model(f.change_ring(K), p)
        cp1 = HyperellipticCurve(f.change_ring(GF(p))).frobenius_polynomial().reverse()
        F1 = hypellfrob(p, max(3, prec), f.lift())
        cp, frob_matrix = from_H1_to_H2(cp1, F1, tensor=kwargs.get('tensor', False))
        frob_matrix = frob_matrix.change_ring(K)
        shift = 0
    elif len(Rf.gens()) == 3 and f.total_degree() == 4 and f.is_homogeneous():
        # Quartic plane curve
        if p < 17:
            prec = max(4, prec)
        else:
            prec = max(3, prec)
        OK = ZpCA(p, prec=prec)
        cp1, F1 = controlledreduction(f,
                                      p,
                                      min_abs_precision=prec,
                                      frob_matrix=True,
                                      threads=1
                                      )
        # change ring to OK truncates precision accordingly
        F1 = F1.change_ring(OK).change_ring(K)
        cp, frob_matrix = from_H1_to_H2(cp1, F1, tensor=kwargs.get('tensor', False))
        shift = 0
    elif len(Rf.gens()) == 4 and f.total_degree() in [4, 5] and f.is_homogeneous():
        shift = 1
        # we will not see the polarization
        # Quartic surface
        if f.total_degree() == 4:
            if p == 3:
                prec = max(5, prec)
            elif p == 5:
                prec = max(4, prec)
            elif p < 43:
                prec = max(3, prec)
            else:
                prec = max(2, prec)
        elif f.total_degree() == 5:
            if p in [3,5]:
                prec = max(7, prec)
            elif p <= 23:
                prec = max(6, prec)
            else:
                prec = max(5, prec)
        OK = ZpCA(p, prec=prec)
        # a rough estimate if it is worth to find a non degenerate mode for f
        if 'find_better_model' in kwargs:
            model = kwargs['find_better_model']
        else:
            model = binomial(3 + (prec - 1)*f.total_degree(), 3) < 6*len(list(f**(prec-1)))
        threads = kwargs.get('threads', ncpus)
        cp, frob_matrix = controlledreduction(f,
                                              p,
                                              min_abs_precision=prec,
                                              frob_matrix=True,
                                              find_better_model=model,
                                              threads=threads
                                              )
        frob_matrix = frob_matrix.change_ring(OK).change_ring(K)
    else:
        raise NotImplementedError("At the moment we only support:\n"
                                  " - Quartic or quintic surfaces\n"
                                  " - Jacobians of quartic curves\n"
                                  " - Jacobians of hyperelliptic curves\n")
    return prec, cp, frob_matrix, shift


def crystalline_obstruction(f, p, precision, over_Qp=False, pedantic=False, **kwargs):
    """
    INPUT:

    - ``f`` -- a polynomial defining the curve or surface.  Note if given an hyperelliptic curve we will try to change the model over Qpbar,
    in order to work with an odd and monic model over Qp.

    - ``p`` -- a prime of good reduction

    - ``precision`` -- a lower bound for the desired precision to run the computations, this increases the time exponentially

    - ``over_Qp`` -- by default False, if True uses the factorization of the cyclotomic polynomials over Qp

    - ``pedantic` -- by default False, if True might inform user of some bound improvements which werent achieved by pure linear algebra arguments

    - ``kwargs`` -- keyword arguments to bypass some computations or to be passed to controlledreduction

    OUTPUT:

    - an upper bound on the number of Tate classes over Qbar

    - a dictionary more information, matching the papers notation, on how one we attained that bound



    EXAMPLES::

    README examples

    Jacobians of hyperelliptic curves with a Weierstrass point over Qp

    Example 5.1

        sage: from crystalline_obstruction import crystalline_obstruction
        sage: f = ZZ['x,y']('x^5 - 2*x^4 + 2*x^3 - 4*x^2 + 3*x - 1 -y^2')
        sage: crystalline_obstruction(f=f, p=31, precision=3) # bounding dim Pic
        (1,
         {'dim Li': [1],
          'dim Ti': [2],
          'factors': [(t - 1, 2)],
          'p': 31,
          'precision': 3,
          'rank T(X_Fpbar)': 2})


    Bounding the geometric dimension of Endomorphism algebra

        sage: f = ZZ['x,y']('x^5 - 2*x^4 + 2*x^3 - 4*x^2 + 3*x - 1 -y^2')
        sage: crystalline_obstruction(f=f, p=31, precision=3, tensor=True) # bounding dim End
        (1,
         {'dim Li': [1],
          'dim Ti': [4],
          'factors': [(t - 1, 4)],
          'p': 31,
          'precision': 3,
          'rank T(X_Fpbar)': 4})

    Example 5.2

        sage: f = ZZ['x,y']('x^5 - 2*x^4 + 7*x^3 - 5*x^2 + 8*x + 3 -y^2')
        sage: crystalline_obstruction(f=f, p=4999, precision=20) # bounding dim Pic
        (2,
         {'dim Li': [2],
          'dim Ti': [2],
          'factors': [(t - 1, 2)],
          'p': 4999,
          'precision': 20,
          'rank T(X_Fpbar)': 2})

    Hyperelliptic curve given in a non-Weierstrass format

        sage: f = ZZ['x,y']('(2*x^6+3*x^5+5*x^4+6*x^3+4*x^2+x) -y*(x^4+x^3+x) -y^2')
        sage: crystalline_obstruction(f=f, p=59, precision=3)
        (3,
         {'dim Li': [1, 2],
          'dim Ti': [3, 6],
          'factors': [(t - 1, 3), (t^2 + t + 1, 3)],
          'p': 59,
          'precision': 3,
          'rank T(X_Fpbar)': 9})

    Jacobians of quartic plane curves
    Example 5.3

        sage: f = ZZ['x,y,z']('x*y^3 + x^3*z - x*y^2*z + x^2*z^2 + y^2*z^2 - y*z^3')
        sage: crystalline_obstruction(f, p=31, precision=3) # bounding dim Pic
        (1,
         {'dim Li': [1],
          'dim Ti': [3],
          'factors': [(t - 1, 3)],
          'p': 31,
          'precision': 3,
          'rank T(X_Fpbar)': 3})

    Product of 3 elliptic curves over x^3 - 3*x - 1

        sage: f=ZZ['x,y,z']('x^3*z + x^2*y*z + x^2*z^2 - x*y^3 - x*y*z^2 - x*z^3 + y^2*z^2')
        sage: crystalline_obstruction(f=f, p=31, precision=3) # bounding dim Pic
        (3,
         {'dim Li': [1, 2],
          'dim Ti': [3, 6],
          'factors': [(t - 1, 3), (t^2 + t + 1, 3)],
          'p': 31,
          'precision': 3,
          'rank T(X_Fpbar)': 9})


    Another gennus 3 plane quartic

        sage: f = ZZ['x,y,z']('x^4+x^2*y^2+2*x^2*y*z-x^2*z^2-6*y^4+16*y^3*z-12*y^2*z^2-16*y*z^3-6*z^4')
        sage: crystalline_obstruction(f=f, p=5003, precision=3) # bounding dim Pic
        (6,
         {'dim Li': [2, 2, 2],
          'dim Ti': [2, 3, 4],
          'factors': [(t + 1, 2), (t - 1, 3), (t^2 + 1, 2)],
          'p': 5003,
          'precision': 3,
          'rank T(X_Fpbar)': 9})
        sage: crystalline_obstruction(f=f, p=5003, precision=3, tensor=True) # bounding dim End
        (9,
         {'dim Li': [2, 3, 4],
          'dim Ti': [4, 6, 8],
          'factors': [(t + 1, 4), (t - 1, 6), (t^2 + 1, 4)],
          'p': 5003,
          'precision': 3,
          'rank T(X_Fpbar)': 18})

    Quartic surfaces

    Example 5.5

        sage: f = ZZ['x,y,z,w']("x^4 + y^4 + z^4 + w^4 + 101^3*x*y*z*w")
        sage: crystalline_obstruction(f, p=101, precision=3)
        (20,
         {'dim Li': [1, 7, 12],
          'dim Ti': [1, 7, 12],
          'factors': [(t - 1, 1), (t - 1, 7), (t + 1, 12)],
          'p': 101,
          'precision': 3,
          'rank T(X_Fpbar)': 20})
        sage: crystalline_obstruction(f=f, p=101, precision=3)
        (19,
         {'dim Li': [1, 6, 12],
          'dim Ti': [1, 7, 12],
          'factors': [(t - 1, 1), (t - 1, 7), (t + 1, 12)],
          'p': 101,
          'precision': 3,
          'rank T(X_Fpbar)': 20})

    Example 5.6

        sage: f = ZZ['x,y,z,w']("y^4 - x^3*z + y*z^3 + z*w^3 + w^4")
        sage: crystalline_obstruction(f=f, p=89, precision=3)
        (4,
         {'dim Li': [1, 0, 3, 0],
          'dim Ti': [1, 1, 4, 4],
          'factors': [(t - 1, 1), (t + 1, 1), (t - 1, 4), (t^4 + 1, 1)],
          'p': 89,
          'precision': 3,
          'rank T(X_Fpbar)': 10})

    Example 5.7

        sage: f = ZZ['x,y,z,w']("x^4 + 2*y^4 + 2*y*z^3 + 3*z^4 - 2*x^3*w- 2*y*w^3")
        sage: crystalline_obstruction(f=f, p=67, precision=3)
        (3,
         {'dim Li': [1, 2],
          'dim Ti': [1, 3],
          'factors': [(t - 1, 1), (t + 1, 3)],
          'p': 67,
          'precision': 3,
          'rank T(X_Fpbar)': 4})

    TESTS::

    Check that precision = 3 is sufficient

        sage: from crystalline_obstruction import crystalline_obstruction
        sage: crystalline_obstruction(f=ZZ['x,y']('y^2-(48*x^8 + 12*x^6 - 22*x^4- 13*x^2 - 2)'),p=107,precision=3)
        (2,
         {'dim Li': [2],
          'dim Ti': [3],
          'factors': [(t - 1, 3)],
          'p': 107,
          'precision': 3,
          'rank T(X_Fpbar)': 3})
        sage: crystalline_obstruction(f=ZZ['x,y']('y^2-(48*x^8 + 12*x^6 - 22*x^4- 13*x^2 - 2)'),p=107,precision=3,tensor=True)
        (2,
         {'dim Li': [2],
          'dim Ti': [6],
          'factors': [(t - 1, 6)],
          'p': 107,
          'precision': 3,
          'rank T(X_Fpbar)': 6})

    """
    if 'cp' in kwargs and 'frob_matrix' in kwargs:
        cp = kwargs['cp']
        frob_matrix = kwargs['frob_matrix']
        shift = kwargs.get('shift', 0)
    else:
        precision, cp, frob_matrix, shift = compute_frob_matrix_and_cp_H2(f, p, precision, **kwargs)
    Rt = PolynomialRing(ZZ, 't')
    t = Rt.gens()[0]
    cp = Rt(cp)
    rank, k, cyc_factor = rank_fieldextension(cp, shift)
    if cyc_factor:
        tate_factor = tate_factor_Zp(cyc_factor.expand())
        max_degree = max(elt.degree() for elt, _ in tate_factor)
        if max_degree > precision - 1:
            warnings.warn('Precision is very likely too low to correctly compute the Tate classes at this prime')
    factor_i, dim_Ti, obsi, obsi_val, dim_Li, _ = upper_bound_tate(cp, frob_matrix, precision, over_Qp=over_Qp, pedantic=pedantic, minors=False)
    res = {}
    res['precision'] = precision
    res['p'] = p
    res['rank T(X_Fpbar)'] = rank
    res['factors'] = []
    if shift > 0:
        res['factors'].append((t - 1, 1))
    # normalize the cyclotomic factors
    for factor, exp in factor_i:
        res['factors'].append((factor(t/p) , exp))
    if over_Qp:
        if shift > 0:
            dim_Li = [[shift]] + dim_Li
            dim_Ti = [[shift]] + dim_Ti
        upper_bound = rank - (sum(sum(dim_Ti, [])) - sum(sum(dim_Li, [])))
    else:
        if shift > 0:
            dim_Li = [shift] + dim_Li
            dim_Ti = [shift] + dim_Ti
        upper_bound = rank - (sum(dim_Ti) - sum(dim_Li))
    res['dim Ti'] = dim_Ti
    res['dim Li'] = dim_Li
    return upper_bound, res






def rank_fieldextension(frob_polynomial, shift=0):
    """
    Return rank, degree of field extension, and the factorization
    of characteristic polynomial, into twisted cyclotomic factors,
    of the Frobenius action on Tate classes factorized

    Input::
        - frob_polynomial, Frobenius polynomial for H^2
        - shift, an integer, 0 by default, accounting for missing cycles,
        for example when frob_polynomial doesn't include the polarization
    """
    p = frob_polynomial.list()[-1].prime_factors()[0]
    rank = shift
    k = 1
    cyc_factorization = []
    for fac, exp in frob_polynomial.factor():
        ki = fac(fac.variables()[0]/p).is_cyclotomic(certificate=True)
        if ki:
            k = lcm(k, ki)
            cyc_factorization.append((fac, exp))
            rank += fac.degree() * exp
    return rank, k, Factorization(cyc_factorization)

def tate_factor_Zp(cyc_factor):
    """
    return the factorization of characteristic polynomial
    of Frobenius on Tate classes in H^2 over Zp
    """
    T = cyc_factor.variables()[0]
    p = cyc_factor.leading_coefficient().prime_factors()[0]
    R = Zp(p)
    factorization = []
    for fac, exp in cyc_factor(T/p).factor():
        for fp, ep in fac.change_ring(R).factor():
            factorization.append((fp(p*T).reverse().monic(), exp*ep))
    res = Factorization(factorization)
    assert res.expand().change_ring(R) == cyc_factor.change_ring(R).reverse().monic(), "%s\n%s" % (res.expand(), cyc_factor.change_ring(R).reverse().monic())
    return res


def upper_bound_tate(cp, frob_matrix, precision, over_Qp=False, pedantic=True, minors=False):
    """
    Return a upper bound for Tate classes over characteristic 0
    TODO: improove documentation
    """
    p = cp.list()[-1].prime_factors()[0]
    # it would be nice to use QpLF
    K = Qp(p, prec=precision+10)
    OK = ZpCA(p, prec=precision)

    # adjust precision
    frob_matrix = matrix(OK, frob_matrix).change_ring(K)

    # get the p-adic eigenvalues
    _, _, cyc_factorization = rank_fieldextension(cp)

    # hacky
    #print(val)
    val = [min(elt.valuation() for elt in col) for col in frob_matrix.columns()]
    projection_cols = frob_matrix.ncols() - val.index(0)
    assert set(val[-projection_cols:]) == {0}
    # P1 = | zero matrix |
    #      | Identity    |
    P1 = zero_matrix(frob_matrix.ncols()-projection_cols, projection_cols).stack(identity_matrix(projection_cols))
    # for the purpose of computing the rank on a projection
    # we want to focus on the first columns
    # thus, we reverse the order of rows and columns
    P1.reverse_rows_and_columns()
    frob_matrix.reverse_rows_and_columns()

    @cached_function
    def frob_power(k):
        if k == 0:
            return identity_matrix(frob_matrix.ncols())
        elif k == 1:
            return frob_matrix
        else:
            return frob_matrix * frob_power(k-1)

    def minors_valuation(M, rank):
        # consider using smith form for minors
        if not minors:
            return None
        if rank == min(M.ncols(), M.nrows()):
            return +Infinity
        elif rank == 0:
            # avoiding bug for 0 matrices
            return min([elt.valuation() for elt in M.list()])
        else:
            return min([elt.valuation() for elt in M.minors(rank + 1)])


    T = matrix(K, 0, frob_matrix.ncols())
    factor_i = []
    dim_Ti = []
    obsi = []
    obsi_val = []
    dim_Li =[]
    dim_Li_val = []
    for cyc_fac, cyc_exp in cyc_factorization:
        factor_i.append((cyc_fac, cyc_exp))
        Ti = matrix(K, 0, frob_matrix.ncols())
        obsij = []
        obsij_val = []
        dim_Tij = []
        dim_Lij = []
        dim_Lij_val = []
        for fac, exp in tate_factor_Zp(cyc_fac):
            # the rows of Tij are a basis for Tij
            # the 'computed' argument avoids echelonizing the kernel basis
            # which might induce some precision loss on the projection
            Tij = fac(frob_matrix).right_kernel_matrix(basis='computed')
            assert Tij.nrows() == fac.degree()*exp*cyc_exp, "We are missing some eigenvalues, increasing the precision should solve this"
            if over_Qp:
                dim_Tij.append(Tij.nrows())
                obs_map = Tij*P1
                rank_obs_ij = obs_map.rank()
                obsij.append(rank_obs_ij)
                obsij_val.append(minors_valuation(obs_map, rank_obs_ij))
                Lijmatrix = matrix(K, Tij.nrows(), 0)
                for ell in range(fac.degree()):
                    Lijmatrix = Lijmatrix.augment(Tij*frob_power(ell).transpose()*P1)
                # Lij = right_kernel(K) subspace of Tij that is invariant under Frob and unobstructed
                Krank = Lijmatrix.rank()
                dim_Lij.append(Tij.nrows() - Krank)
                dim_Lij_val.append(minors_valuation(Lijmatrix, Krank))
                assert dim_Lij[-1] % fac.degree() == 0

            Ti = Ti.stack(Tij)
        T = T.stack(Ti)
        if over_Qp:
            dim_Ti.append(dim_Tij)
            obsi.append(obsij)
            obsi_val.append(obsij_val)
            dim_Li.append(dim_Lij)
            dim_Li_val.append(dim_Lij_val)
        else:
            obs_map = Ti*P1
            assert Ti.nrows() == cyc_fac.degree()*cyc_exp
            dim_Ti.append(Ti.nrows())
            rank_obs_i = obs_map.rank()
            obsi.append(rank_obs_i)
            obsi_val.append(minors_valuation(obs_map, rank_obs_i))
            Limatrix = matrix(K, Ti.nrows(), 0)
            for ell in range(0,cyc_fac.degree()):
                Limatrix = Limatrix.augment(Ti*frob_power(ell).transpose()*P1)
            # Li = right_kernel(K) subspace of Tij that is invariant under Frob and unobstructed
            Krank = Limatrix.rank()
            dim_Li.append(Ti.nrows() - Krank)
            dim_Li_val.append(minors_valuation(Limatrix, Krank))
            if dim_Li[-1] % cyc_fac.degree()  != 0:
                old_dim = dim_Li[-1]
                deg = cyc_fac.degree()
                new_dim = dim_Li[-1] = deg * (old_dim // deg)
                if pedantic:
                    warnings.warn("rounding dimension of Li from %d to %d for cyc_factor = %s" % (old_dim, new_dim, cyc_fac))
    assert T.nrows() == sum(fac.degree()*exp for fac, exp in cyc_factorization)

    return factor_i, dim_Ti, obsi, obsi_val, dim_Li, dim_Li_val
