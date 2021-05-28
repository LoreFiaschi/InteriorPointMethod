include("../src/ipqp.jl");
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl");
import LinearAlgebra
using .BAN


# DISJUNCT PROBLEM FROM SLIDE 26 PACK 3

Q = [-I zeros(Ban, 2,4); zeros(Ban, 4,6)];

A = [-1  1  1  0  0  0;
     -1  1  0 -1  0  0;
      1  1  0  0 -1  0;
      1  1  0  0  0  1;
	];

A = convert(Matrix{Ban}, A)

b = [ 1, -1, 3, 5];

b = convert(Vector{Ban},b)

c = [ 2, 1, 0, 0, 0, 0].*η;


tol=1e-8;
verbose = false;
genLatex = false;

sol = ipqp(A,b,c,Q, tol; maxit=20, verbose=verbose, genLatex=genLatex, slack_var=3:6);
nothing

