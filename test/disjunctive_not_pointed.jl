include("../src_debug/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

#
# x1 => 2 && x2 <= -4
#
# -x1 + x2 <= 1
#

M = α;

#A = [I -I -I zeros(Ban,4,2);
#	zeros(Ban,1,12) ones(Ban,1,2);
#	zeros(Ban,8,4) I M.*[-ones(Ban,4,1) zeros(Ban,4,1); zeros(Ban,4,1) -ones(Ban,4,1)]];
#
#A_dom = [zeros(Ban,3,4) [1 -1 zeros(Ban,1,6) -2 0;
#						 0  0 1 -1 zeros(Ban,1,4) 4 0;
#						zeros(Ban,1,4) 1 -1 -1 1 0 1]];
#
#b = [zeros(Ban,4); 1; zeros(Ban, 11)];
#
#t = [zeros(Int64,5); -ones(Int64,8); 1; -1; 1];
#
##ge_idx = findall(x->x>0, t);
#A[ge_idx, :] .*= -1;
#b[ge_idx] .*= -1;
#
#A = [A [zeros(Ban,5,size(A,1)-5); I]]
#

A = [zeros(Ban,1,8) ones(Ban,1,2);
	 I M.*[-ones(Ban,4,1) zeros(Ban,4,1); zeros(Ban,4,1) -ones(Ban,4,1)]];

A_dom = [-1  1  0  0 zeros(Ban, 1, 4) 2  0;
		  0  0  1 -1 zeros(Ban, 1, 4) 4  0;
		  zeros(Ban, 1, 4) -1 1 1 -1  0 -1];

A = [A; A_dom];

A = [A [zeros(Ban,1,size(A,1)-1); I]]

b = [1; zeros(Ban, size(A,1)-1)];

c = [1+η;  -1+η;  -1+η;  1+η; zeros(Ban,size(A,2)-4)];
Q = zeros(Ban, size(A,2), size(A,2));

tol = 1e-8;
genLatex = false;
verbose = false;

sol = ipqp(A,b,c,Q, tol; maxit=10, verbose=verbose, genLatex=genLatex, slack_var=5:size(A,2));# slack_var=11:21);

print("\tSolution: "); 
println([sol.x[1]-sol.x[2]; sol.x[3]-sol.x[4]]);
print("\tDisjoint flag: "); println(sol.x[end-1:end]);
