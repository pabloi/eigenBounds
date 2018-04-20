# eigenBounds
Code implementing the computation of bounds on eigenvalues/singular values of a matrix.
Work in progress.


# TO DO on Gao implementation:
1) Check if PRbound is implemented properly
2) See if periodic/pseudo-periodic signals can also use the factorization of covariance trick to get a bound
(periodic signals naturaly have a dimensionality bound, but they do not fit within the given framework)
3) See if results can be extended to  circulant matrices 
(i.e. periodic auto-correlations which do not taper off as time goes to infinity)

# TO DO on psd matrix bound
1) Given a lower bound on sum of largest eigenvalues and rank, create the best possible bound -> DONE
2) Given a lower bound on sum of largest, and an upper bound on sum of lowest, create best possible bound: i.e. take max(lowerBound,tr(C)-upperBound), and if result is not a concave bound (it need not be), then generate the minimum cocave curve that is above the max() curve.


