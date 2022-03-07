include("../src_debug/ipqp_u.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

#
# x1 => 2 && x2 <= -4
#
# -x1 + x2 <= 1
#
#
# transformed into
#
#
#
#

M = α;

A = [zeros(Ban,1,4) one(Ban);
	 Matrix{Ban}(I,4,4) M.*[ones(Ban,2,1); -ones(Ban,2,1)];
	-Matrix{Ban}(I,2,2) zeros(Ban, 2, 2) M.*ones(Ban,2,1)];

A_dom = [zeros(Ban, 1, 2) -1  0  2;
		 zeros(Ban, 1, 2)  0 -1  4;
		 -1  1 zeros(Ban, 1, 2)  1];

A = [A; A_dom];

A = [A I]

b = [1; M; M; zero(Ban); zero(Ban); M; M; zero(Ban); zero(Ban); 1];

c = [1;  -1; 1+η;  1+η; zeros(Ban,size(A,2)-4)];
#Q = zeros(Ban, size(A,2), size(A,2));
Q = [I.*η zeros(Ban, 2, size(A,2)-2); zeros(Ban, size(A,2)-2, size(A,2))]

tol = 1e-8;
genLatex = false;
verbose = false;

sol = ipqp(A,b,c,Q, tol; maxit=25, verbose=verbose, genLatex=genLatex, slack_var=6:size(A,2), bounded_variables=3:size(A,2));# slack_var=11:21);

print("\tSolution: "); 
println([sol.x[1], sol.x[2], sol.x[3], sol.x[4]]);
print("\tDisjoint flag: "); println(sol.x[5]);
