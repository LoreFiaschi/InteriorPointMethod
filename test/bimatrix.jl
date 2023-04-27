#include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s2_isbits.jl")
include("../src/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN
using LinearAlgebra

# function to convert a bimatrix game into a QP (assume A and B represents the players losses)
function bm2qp(A::Matrix, B::Matrix)
	
	dim_x, dim_y = size(A);
	dim_xy = dim_x + dim_y;
	M = [zeros(dim_x, dim_x) A; B' zeros(dim_y, dim_y)];
	R = [M -I];
	Q = 2*[M zeros(dim_xy, dim_xy); zeros(dim_xy, 2*dim_xy)];
	c = [-ones(dim_xy); zeros(dim_xy)];

	return Q, c, R;
end


tol=1e-8;
genLatex = true;
verbose = false;

A = [1 1; 0 3];
B = [3 0; 1 2];

#A = [1 0 0; 0 0 3; 0 2 0];
#B = [1 0 0; 0 2 0; 0 0 3];

#sol = [1.66666; 2; 2.16666; 1.83333; 0.5; 0.5; 0.5; 0.5]


Q,c,R = bm2qp(A,B);

Q = convert(Matrix{Ban}, Q);
c = convert(Vector{Ban}, c);
R = convert(Matrix{Ban}, R);

sol = ipqp(R,ones(Ban, size(R,1)),c,Q, tol; maxit=20, verbose=verbose, genLatex=genLatex, slack_var=5:8);

nothing

