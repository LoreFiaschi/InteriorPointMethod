include("../src/ipqp.jl");
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl");
import LinearAlgebra
using .BAN


# DISJUNCT PROBLEM FROM SLIDE 26 PACK 3

Q = [ 0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0  -α^2 0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
      0  0   0   0  0  0  0  0  0  0  0  0;
    ];

A = [ 1  0  α  1  0  0  0  0  0  0  0  0;
      1  0 -α  0 -1  0  0  0  0  0  0  0;
      0  1  α  0  0  1  0  0  0  0  0  0;
      0  1 -α  0  0  0 -1  0  0  0  0  0;
      1  0 -α  0  0  0  0  1  0  0  0  0;
      1  0  α  0  0  0  0  0 -1  0  0  0;
      0  1 -α  0  0  0  0  0  0  1  0  0;
      0  1  α  0  0  0  0  0  0  0 -1  0;
      0  0  1  0  0  0  0  0  0  0  0  1;
    ];

A = convert(Matrix{Ban}, A);

b = [ α+3, -α, α+4, -α, 9, 5, 6, 4, 1 ];

# maximize -1 -1
#c = [ 1, 1, -2α, -2α, 0, 0, 0, 0, 0, 0, 0, 0 ];

# maximize 1 1
c = [ -1, -1, α^2, 0, 0, 0, 0, 0, 0, 0, 0, 0];


tol=1e-8;
verbose = false;
genLatex = false;

sol = ipqp(A,b,c,Q, tol; maxit=15, verbose=verbose, genLatex=genLatex, slack_var=5:12);
nothing

