# crystalline_obstruction

This package computes an approximation the crystalline obstruction map on its space of Tate classes.
This gives rigorous upper bounds on the middle geometric Picard number of a given projective hypersurface or the geometric Picard number of Jacobian variety.

For more details see:

- "Explicit computation of the obstruction to lifting Tate classes from positive characteristic" by Edgar Costa and Emre Can Sertoz.


To compute a p-adic approximation for the Frobenius matrix we rely on the libraries [controlledreduction](https://github.com/edgarcosta/controlledreduction), and ([hypellfrob](https://web.maths.unsw.edu.au/~davidharvey/code/hypellfrob/)).
At the moment we support:

- Jacobians hyperelliptic curves with a Weierstrass point over Qp ([hypellfrob](https://web.maths.unsw.edu.au/~davidharvey/code/hypellfrob/))

- Jacobians of quartic curves ([controlledreduction](https://github.com/edgarcosta/controlledreduction))

- quartic and quintic surfaces ([controlledreduction](https://github.com/edgarcosta/controlledreduction))




# Installing

```
sage -pip install --upgrade git+https://github.com/edgarcosta/crystalline_obstruction.git@master#egg=crystalline_obstruction
```
If you don't have permissions to install it system wide, please add the flag ``--user`` to install it just for you.
```
sage -pip install --user --upgrade git+https://github.com/edgarcosta/crystalline_obstruction.git@master#egg=crystalline_obstruction
```

# Examples

## Jacobians of hyperelliptic curves with a Weierstrass point over Qp

```
sage: from crystalline_obstruction import crystalline_obstruction
sage: #Example 5.1
....: crystalline_obstruction(f=ZZ['x,y']('x^5 - 2*x^4 + 2*x^3 - 4*x^2 + 3*x - 1 -y^2'),
....:                         p=31,
....:                         precision=3)
(1,
 {'precision': 3,
  'p': 31,
  'rank T(X_Fpbar)': 2,
  'factors': [(t - 1, 2)],
  'dim Ti': [2],
  'dim Li': [1]})
sage: #Example 5.2, with prec = 100 takes about 3 minutes
sage: crystalline_obstruction(f=ZZ['x,y']('x^5 - 2*x^4 + 7*x^3 - 5*x^2 + 8*x + 3 -y^2'),
....:                         p=4999,
....:                         precision=20)
(2,
 {'precision': 20,
  'p': 4999,
  'rank T(X_Fpbar)': 2,
  'factors': [(t - 1, 2)],
  'dim Ti': [2],
  'dim Li': [2]})
sage: crystalline_obstruction(f=ZZ['x,y']('(2*x^6+3*x^5+5*x^4+6*x^3+4*x^2+x) -y*(x^4+x^3+x) -y^2'),
....:                         p=59,
....:                         precision=4)
(3,
 {'precision': 4,
  'p': 59,
  'rank T(X_Fpbar)': 9,
  'factors': [(t - 1, 3), (t^2 + t + 1, 3)],
  'dim Ti': [3, 6],
  'dim Li': [1, 2]})
```

## Quartic Plane curves
```
sage: #Example 5.3
sage: crystalline_obstruction(f=ZZ['x,y,z']('x*y^3 + x^3*z - x*y^2*z + x^2*z^2 + y^2*z^2 - y*z^3'),
....:                         p=31,
....:                         precision=3)
(1,
 {'precision': 3,
  'p': 31,
  'rank T(X_Fpbar)': 3,
  'factors': [(t - 1, 3)],
  'dim Ti': [3],
  'dim Li': [1]})
sage: # Product of 3 elliptic curves over x^3 - 3*x - 1
sage: crystalline_obstruction(f=ZZ['x,y,z']('x^3*z + x^2*y*z + x^2*z^2 - x*y^3 - x*y*z^2 - x*z^3 + y^2*z^2'),
....:                         p=31,
....:                         precision=5)
(3,
 {'precision': 5,
  'p': 31,
  'rank T(X_Fpbar)': 9,
  'factors': [(t - 1, 3), (t^2 + t + 1, 3)],
  'dim Ti': [3, 6],
  'dim Li': [1, 2]})
```

## Quartic surfaces
```
sage: # Example 5.5
sage: crystalline_obstruction(ZZ['x,y,z,w']("x^4 + y^4 + z^4 + w^4 + 101^3*x*y*z*w"),
....:                         p=101,
....:                         precision=3)
(20,
 {'precision': 3,
  'p': 101,
  'rank T(X_Fpbar)': 20,
  'factors': [(t - 1, 1), (t - 1, 7), (t + 1, 12)],
  'dim Ti': [1, 7, 12],
  'dim Li': [1, 7, 12]})
sage: # Example 5.5
sage: crystalline_obstruction(ZZ['x,y,z,w']("x^4 + y^4 + z^4 + w^4 + 101^3*x*y*z*w"),
....:                         p=101,
....:                         precision=4)
(19,
 {'precision': 4,
  'p': 101,
  'rank T(X_Fpbar)': 20,
  'factors': [(t - 1, 1), (t - 1, 7), (t + 1, 12)],
  'dim Ti': [1, 7, 12],
  'dim Li': [1, 6, 12]})
sage: # Example 5.6
sage: crystalline_obstruction(ZZ['x,y,z,w']("y^4 - x^3*z + y*z^3 + z*w^3 + w^4"),
....:                         p=89,
....:                         precision=3)
(4,
 {'precision': 3,
  'p': 89,
  'rank T(X_Fpbar)': 10,
  'factors': [(t - 1, 1), (t + 1, 1), (t - 1, 4), (t^4 + 1, 1)],
  'dim Ti': [1, 1, 4, 4],
  'dim Li': [1, 0, 3, 0]})
sage: # Example 5.7
sage: crystalline_obstruction(ZZ['x,y,z,w']("x^4 + 2*y^4 + 2*y*z^3 + 3*z^4 - 2*x^3*w- 2*y*w^3"),
....:                         p=67,
....:                         precision=3)
(3,
 {'precision': 3,
  'p': 67,
  'rank T(X_Fpbar)': 4,
  'factors': [(t - 1, 1), (t + 1, 3)],
  'dim Ti': [1, 3],
  'dim Li': [1, 2]})
```



# Citing this code

Please cite the following preprint if this code has been helpful in your research:

- "Explicit computation of the obstruction to lifting Tate classes from positive characteristic" by Edgar Costa and Emre Can Sertoz.

Preprint available at arXiv:???
