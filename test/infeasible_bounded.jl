include("../src_debug/ipqp_u.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
#include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s3_isbits.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

#
# min - x1 - x2
#
# x1 >= 1 && x2 >= 1 && x1 <= 2 && x2 <= 2
#
# x1 - x2 >= 1 && x1 - x2 <= -1 (unbounded vars)

# vars: x12 x22 x11 x12 y

M = α;

A = [zeros(Ban,1,4) one(Ban);
	 Matrix{Ban}(I,4,4) M.*[ones(Ban,2,1); -ones(Ban,2,1)];
	-Matrix{Ban}(I,2,2) zeros(Ban, 2, 2) M.*ones(Ban,2,1)];

A_dom = [-1  1 zeros(Ban, 1, 2) -1;
		  1 -1 zeros(Ban, 1, 2) -1;
		  zeros(Ban, 2, 2)  Matrix{Ban}(I,2,2)    [-2; -2];
		  zeros(Ban, 2, 2) -Matrix{Ban}(I,2,2)    [ 1;  1] ];

A = [A; A_dom];
A = [A I]

b = [1; M; M; zero(Ban); zero(Ban); M; M; -1; -1; 0; 0; 0; 0];

c = [-1; -1; -1; -1; zeros(Ban,size(A,2)-4)];

Q = zeros(Ban, size(A,2), size(A,2));

tol = 1e-8;
genLatex = false;
verbose = false;

sol = ipqp(A,b,c,Q, tol; maxit=16, verbose=verbose, genLatex=genLatex, slack_var=6:size(A,2), bounded_variables=3:size(A,2));

nothing

#print("\tSolution: "); 
#println([x[1], x[2], x[3], x[4]]);
#print("\tDisjoint flag: "); println(x[5]);
