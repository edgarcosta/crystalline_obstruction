from sage.all import ZZ # cannot seem to import this from the right place
from sage.functions.other import binomial
from sage.rings.padics.factory import Qp, Zp, ZpCA
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.arith.functions import lcm
from sage.matrix.special import zero_matrix, identity_matrix, matrix, companion_matrix
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.structure.factorization import Factorization
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from pycontrolledreduction import controlledreduction
import multiprocessing
ncpus = multiprocessing.cpu_count()
import warnings
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob
from sage.misc.cachefunc import cached_function


def from_H1_to_H2(cp1, F1):
    """
    Input:
        cp1 - characteristic polynomial of Frob acting on H^1
        F1 - Frob approximation of the action H^1
    Output:
        cp2 - characteristic polynomial of Frob acting on H^2
        F2 - Frob approximation of the action H^2

    """
    # deal with Frob approximation
    wedge_basis = [(i,j) for i in range(F1.nrows()) for j in range(i+1, F1.nrows())]
    M = FiniteRankFreeModule(F1.base_ring(), F1.nrows())
    M.basis('e')
    F1e = [M(elt) for elt in F1.columns()]
    F2e = [F1e[i].wedge(F1e[j]) for i, j in wedge_basis]
    F2 = matrix([[elt.comp()[i,j] for elt in F2e] for i,j in wedge_basis])
    # deduce the correct characteristic polynomial
    assert cp1[0] == 1, "The constant term of the characteristic polynomial must be 1"
    cprev = cp1.reverse()
    F1_cp = companion_matrix(cprev)
    assert F1.nrows() == F1_cp.nrows()
    M = FiniteRankFreeModule(F1_cp.base_ring(), F1.nrows())
    M.basis('e')
    F1e_cp = [M(elt) for elt in F1_cp.columns()]
    F2e_cp = [F1e_cp[i].wedge(F1e_cp[j]) for i, j in wedge_basis]
    F2_cp = matrix([[elt.comp()[i,j] for elt in F2e_cp] for i,j in wedge_basis])
    cp2 = F2_cp.charpoly(var = cp1.variables()[0]).reverse()
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


def compute_frob_matrix_and_cp_H2(f, p, prec):
    """
    Return a p-adic matrix approximating the action of Frob on H^2 of a surface or abelian variety,
    and its characteristic polynomial over ZZ
    Input:
        - f, or fs defining the curve or surface
        - p, prime
        - prec, a lower bound for the desired precision to run the computations, this increases the time exponentially

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
        F1 = hypellfrob(p, prec, f.lift())
        cp, frob_matrix = from_H1_to_H2(cp1, F1)
        frob_matrix = frob_matrix.change_ring(K)
        shift = 0
    elif len(Rf.gens()) == 3 and f.total_degree() == 4 and f.is_homogeneous():
        # Quartic plane curve
        if p < 17:
            prec = max(3, prec)
        else:
            prec = max(2, prec)
        OK = ZpCA(p, prec=prec)
        cp1, F1 = controlledreduction(f,
                                      p,
                                      min_abs_precision=prec,
                                      frob_matrix=True,
                                      threads=ncpus
                                      )
        # change ring to OK truncates precision accordingly
        F1 = F1.change_ring(OK).change_ring(K)
        cp, frob_matrix = from_H1_to_H2(cp1, F1)
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
        model = binomial(3 + (prec - 1)*f.total_degree(), 3) < 6*len(list(f**(prec-1)))
        cp, frob_matrix = controlledreduction(f,
                                              p,
                                              min_abs_precision=prec,
                                              frob_matrix=True,
                                              find_better_model=model,
                                              threads=ncpus
                                              )
        frob_matrix = frob_matrix.change_ring(OK).change_ring(K)
    else:
        raise NotImplementedError("At the moment we only support:\n"
                                  " - Quartic or quintic surfaces\n"
                                  " - Jacobians of quartic curves\n"
                                  " - Jacobians of hyperelliptic curves\n")
    return prec, cp, frob_matrix, shift


def crystalline_obstruction(f, p, precision, over_Qp=False, **kwargs):
    """
    Input:
        - f, or fs defining the curve or surface
        - p, prime
        - prec, a lower bound for the desired precision to run the computations, this increases the time exponentially

    Output:
        - `prec`, the minimum digits absolute precision for approximation of the Frobenius
        - the number of Tate classes
        - a lower bound on the rank of the obstruction map

    Note: if given a univariate polynomial (or a pair), we will try to change the model over Qpbar,
    in order to work with an odd and monic model
    """
    if 'cp' in kwargs and 'frob_matrix' in kwargs:
        cp = kwargs['cp']
        frob_matrix = kwargs['frob_matrix']
        shift = kwargs.get('shift', 0)
    else:
        precision, cp, frob_matrix, shift = compute_frob_matrix_and_cp_H2(f, p, precision)
    Rt = PolynomialRing(ZZ, 't')
    t = Rt.gens()[0]
    cp = Rt(cp)
    rank, k, cyc_factor = rank_fieldextension(cp, shift)
    tate_factor = tate_factor_Zp(cyc_factor.expand())
    max_degree = max(elt.degree() for elt, _ in tate_factor)
    if max_degree > precision - 1:
        warnings.warn('Precision is very likely too low to correctly compute the Tate classes at this prime')
    factor_i, dim_Ti, obsi, obsi_val, dim_Li, dim_Li_val = upper_bound_tate(cp, frob_matrix, precision, over_Qp=over_Qp)
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


def upper_bound_tate(cp, frob_matrix, precision, over_Qp=False):
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
            Tij = fac(frob_matrix).right_kernel_matrix()
            assert Tij.nrows() == fac.degree()*exp*cyc_exp
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
            assert dim_Li[-1] % cyc_fac.degree()  == 0
    assert T.nrows() == sum(fac.degree()*exp for fac, exp in cyc_factorization)

    return factor_i, dim_Ti, obsi, obsi_val, dim_Li, dim_Li_val
