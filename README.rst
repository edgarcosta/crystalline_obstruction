=============
crystalline_obstruction
=============

This package computes an approximation the crystalline obstruction map on its space of Tate classes.
This gives rigorous upper bounds on the middle geometric Picard number of a given projective hypersurface or the geometric Picard number of Jacobian variety.

For more details see

  - Explicit computation of the obstruction to lifting Tate classes from positive characteristic (under preparation) by Edgar Costa and Emre Can Sertoz.

To compute a p-adic approximation for the Frobenius matrix we rely on the libraries controlledreduction_, and hypellfrob_.
At the moment we support:
  - quartic or quintic surfaces (controlledreduction_)
  - Jacobians of quartic curves (controlledreduction_)
  - Jacobians hyperelliptic curves with a Weierstrass point over Qp (hypellfrob_)

.. _controlledreduction: https://github.com/edgarcosta/controlledreduction
.. _hypellfrob: https://web.maths.unsw.edu.au/~davidharvey/code/hypellfrob/

============
Install
============

::

  sage -pip install --upgrade git+https://github.com/edgarcosta/crystalline_obstruction.git@master#egg=crystalline_obstruction


If you don't have permissions to install it system wide, please add the flag ``--user`` to install it just for you.

::

  sage -pip install --user --upgrade git+https://github.com/edgarcosta/crystalline_obstruction.git@master#egg=crystalline_obstruction


============
Examples
============

::

    sage: from crystalline_obstruction import improve_bound_rank
    sage: improve_bound_rank(f=ZZ['x,y,z,w']('y^4-x^3*z+y*z^3+z*w^3+w^4'),
    ....:                               p=31,
    ....:                               precision=3)

    {'precision': 3, 'p': 31, 'rank T(X_Fpbar)': 4, 'rank obs': 0}
    sage:
    sage: improve_bound_rank(f=ZZ['x,y,z,w']('y^4-x^3*z+y*z^3+z*w^3+w^4'),
    ....:                               p=89,
    ....:                               precision=3,
    ....:                               per_cyclic=True)
    {'precision': 3,
     'p': 89,
     'rank T(X_Fpbar)': 10,
     'rank obs|Ti': [1, 1, 1, 1, 1, 1],
     'dim Ti': [1, 1, 1, 1, 1, 4],
     'sum rank  obs|Ti': 6,
     'rank obs': 1}
    sage: improve_bound_rank(f=ZZ['x,y,z']('-x*y^3 + x^3*z + x^2*y*z + x^2*z^2 - x*y*z^2 + y^2*z^2 - x*z^3'),
    ....:                               p=31,
    ....:                               precision=3)
    {'precision': 3, 'p': 31, 'rank T(X_Fpbar)': 9, 'rank obs': 3}
    sage: improve_bound_rank(f=ZZ['x,y,z']('-x*y^3 + x^3*z + x^2*y*z + x^2*z^2 - x*y*z^2 + y^2*z^2 - x*z^3'),
    ....:                               p=31,
    ....:                               precision=3,
    ....:                               per_cyclic=True)
    {'precision': 3,
     'p': 31,
     'rank T(X_Fpbar)': 9,
     'rank obs|Ti': [2, 2, 2],
     'dim Ti': [3, 3, 3],
     'sum rank  obs|Ti': 6,
     'rank obs': 3}
    sage: improve_bound_rank(f=[ZZ['x']('2*x^6+3*x^5+5*x^4+6*x^3+4*x^2+x'),
    ....:                                  ZZ['x']('x^4+x^3+x')],
    ....:                               p=59,
    ....:                               precision=3,
    ....:                               per_cyclic=True)
    {'precision': 3,
     'p': 59,
     'rank T(X_Fpbar)': 9,
     'rank obs|Ti': [2, 3],
     'dim Ti': [3, 6],
   'sum rank  obs|Ti': 5,
   'rank obs': 3}
    sage: #Example 5.1
    ....: bound_rank.improve_bound_rank(f=ZZ['x']('x^5 - 2*x^4 + 2*x^3 - 4*x^2 + 3*x - 1'),
    ....:                               p=31,
    ....:                               precision=3,
    ....:                               per_cyclic=True)
    {'precision': 3,
     'p': 31,
     'rank T(X_Fpbar)': 2,
     'rank obs|Ti': [1],
     'dim Ti': [2],
     'sum rank  obs|Ti': 1,
     'rank obs': 1}
    sage: #Example 5.2, prec should be 100, but takes some time
    ....: bound_rank.improve_bound_rank(f=ZZ['x']('x^5 - 2*x^4 + 7*x^3 - 5*x^2 + 8*x + 3'),
    ....:                               p=4999,
    ....:                               precision=20,
    ....:                               per_cyclic=True)
    {'precision': 20,
     'p': 4999,
     'rank T(X_Fpbar)': 2,
     'rank obs|Ti': [0],
     'dim Ti': [2],
     'sum rank  obs|Ti': 0,
     'rank obs': 0}
    sage: #Example 5.3
    ....: bound_rank.improve_bound_rank(f=ZZ['x,y,z']('x*y^3 + x^3*z - x*y^2*z + x^2*z^2 + y^2*z^2 - y*z^3'),
    ....:                               p=31,
    ....:                               precision=3,
    ....:                               per_cyclic=True)
    {'precision': 3,
     'p': 31,
     'rank T(X_Fpbar)': 3,
     'rank obs|Ti': [2],
     'dim Ti': [3],
     'sum rank  obs|Ti': 2,
     'rank obs': 2}

    {'precision': 3, 'p': 31, 'rank T(X_Fpbar)': 4, 'rank obs': 0}

============
Citing this code
============

Please cite the following preprint if this code has been helpful in your research:

???

Preprint available at arXiv:???
