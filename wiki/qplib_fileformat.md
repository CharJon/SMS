# Fileformat of QPlib

See: https://static-content.springer.com/esm/art%3A10.1007%2Fs12532-018-0147-4/MediaObjects/12532_2018_147_MOESM1_ESM.pdf
(Appendix B)

In general, the unconstrained version of Quadratic Programs looks like:

```latex
min/max $1/2 x^T * Q^0 * x + b^0 * x + q^0$
$x \in Z$
```

Naming of filetypes works like "OVC", where:

- O: type of objective function considered (Q is quadratic)
- V: types of variables (B is binary)
- C: type of constraints (N is none)

For now, we only care about the **QBN**, which is the same as **QBB**:

```
name:           problem name (character string)
type:           problem type (character string) 
sense:          (minimize | maximize)
n:              number of variables (integer)
n^Q^0:          number of nonzeros (integer) in lower triangle of Q^0
h k Q^0_{hk}:   row and column indices (integers) and value (real) for each
b^0_d:          default value (real) for entries in b^0
n^b^0:          number of non-default entries (integer) in b^0
j b_j^0:        index (integer) and value (real) for each non-default term in b^0, if n^b^0 > 0, one pair per line
q^0:            constant part of objective function
c_inf:          value (real) for infinity for constraint or variable bounds—any, bound greater than or equal to this in, absolute value, is infinite
x^0_d:          default value (real) for the components of the starting point x^0 for the variables x
n^x^0:          number of non-default starting entries (integer) in x
i x_i^0:        index (integer) and value (real) for each non-default starting point x^0, if n^x^0 > 0, one pair per line
```
