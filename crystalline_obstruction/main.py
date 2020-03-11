from sage.functions.other import binomial
from sage.rings.padics.factory import Qp, Zp, ZpCA
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import Infinity
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
            raise NotImplementedError("HyperellipticCurve doesn't have a weirstrass point over Qp")
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
    if isinstance(f, tuple) or isinstance(f, list):
        # assume is hyperelliptic
        # and convert to weirstrass eqn
        fn, hn = f
        f = 4*fn + hn**2
    Rf = f.parent()
    if len(Rf.gens()) == 1:
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


def crystalline_obstruction(f, p, precision, per_cyclic=False):
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
    precision, cp, frob_matrix, shift = compute_frob_matrix_and_cp_H2(f, p, precision)
    rank, k, cyc_factor = rank_fieldextension(cp, shift)
    tate_factor = tate_factor_Zp(cyc_factor)
    max_degree = max(elt.degree() for elt, _ in tate_factor)
    if max_degree > precision - 1:
        warnings.warn('Precision is very likely too low to correctly compute the Tate classes at this prime')
    bounds = upper_bound_tate(cp, frob_matrix, precision, per_cyclic=per_cyclic)
    res = {}
    res['precision'] = precision
    res['p'] = p
    res['rank T(X_Fpbar)'] = rank
    if per_cyclic:
        b, b_val, b_local, b_val_local = bounds
        res['rank obs|Ti'] = b_local
        res['dim Ti'] = [fac.degree()*exp for fac, exp in tate_factor]
        res['sum rank  obs|Ti'] = sum(b_local)
        upper_bound = rank - sum(b_local)
    else:
        b, b_val = bounds
        upper_bound = rank - b
    res['rank obs'] = b
    return upper_bound, res






def rank_fieldextension(frob_polynomial, shift=0):
    """
    Return rank, degree of field extension, and the characteristic polynomial
    of Frobenius on Tate classes.

    Input::
        - frob_polynomial, Frobenius polynomial for H^2
        - shift, an integer, 0 by default, accounting for missing cycles,
        for example when frob_polynomial doesn't include the polarization
    """
    p = frob_polynomial.list()[-1].prime_factors()[0]
    rank = shift
    k = 1
    cyc_factor = 1
    for fac, exp in frob_polynomial.factor():
        ki = fac(fac.variables()[0]/p).is_cyclotomic(certificate=True)
        if ki:
            k = lcm(k, ki)
            cyc_factor *= fac**exp
            rank += fac.degree() * exp
    return rank, k, cyc_factor

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


def upper_bound_tate(cp, frob_matrix, precision, per_cyclic=False):
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
    _, _, cyc_factor = rank_fieldextension(cp)
    tate_factor = tate_factor_Zp(cyc_factor)
    dim_Tn = sum(fac.degree()*exp for fac, exp in tate_factor)

    # hacky
    #print(val)
    val = [min(elt.valuation() for elt in col) for col in frob_matrix.columns()]
    projection_cols = frob_matrix.ncols() - val.index(0)
    assert set(val[-projection_cols:]) == {0}
    P1 = zero_matrix(frob_matrix.ncols()-projection_cols, projection_cols).stack(identity_matrix(projection_cols))
    # for the purpose of computing the rank on a projection
    # we want to focus on the first columns
    # thus, we reverse the order of rows and columns
    P1.reverse_rows_and_columns()
    frob_matrix.reverse_rows_and_columns()

    M = matrix(K, 0, frob_matrix.ncols())
    b_local = []
    bval_local = []
    for fac, exp in tate_factor:
        #print(fac.change_ring(R), exp)
        Mlocal = fac(frob_matrix)
        kerlocal = Mlocal.right_kernel_matrix()
        assert kerlocal.nrows() == fac.degree()*exp
        M = M.stack(kerlocal)
        if per_cyclic:
            obs = kerlocal*P1
            b = obs.rank()
            if b == min(obs.ncols(), obs.nrows()):
                b_val = +Infinity
            elif b == 0:
                # avoiding bug for 0 matrices
                b_val = min([elt.valuation() for elt in obs.list()])
            else:
                b_val = min([elt.valuation() for elt in obs.minors(b + 1)])
            b_local.append(b)
            bval_local.append(b_val)
        #print(kerlocal)
    assert dim_Tn == M.nrows()


    obs = M*P1
    b = obs.rank()
    #print(obs.transpose())
    if b == min(obs.ncols(), obs.nrows()):
        b_val = +Infinity
    elif b == 0:
        # avoiding bug for 0 matrices
        b_val = min([elt.valuation() for elt in obs.list()])
    else:
        b_val = min([elt.valuation() for elt in obs.minors(b + 1)])
    if per_cyclic:
        return b, b_val, b_local, bval_local
    else:
        return b, b_val
