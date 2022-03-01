include("../src_debug/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

# max -x1 + x_2
#
# -x1+ x2 <=  0    y=0
# -x1+ x2 <=  0    y=1
# -x1-2x2 <= -2

function load_param(experiment)

	# standard
	if experiment == 0
		return [-1;  1;  1; -1; -1;  1;  1; -1; zeros(Ban,2+11)], zeros(Ban,21,21);
	# min l1 norm
	else
		return [1+η;  -1+η;  -1+η;  1+η;  1+η;  -1+η;  -1+η;  1+η; zeros(Ban,2+11)], zeros(Ban,21,21);
	end
end

experiment = 1;
M = α;

#A = [I -I -I zeros(Ban,4,2);
#     zeros(Ban,1,12) ones(Ban,1,2);
#	 zeros(Ban,8,4) I M.*[-ones(Ban,4,1) zeros(Ban,4,1); zeros(Ban,4,1) -ones(Ban,4,1)]];
#
#A_dom = [zeros(Ban,3,4) [-1  1  1 -1 zeros(Ban, 1, 4)  0  0;
#                         zeros(Ban, 1, 4) -1  1  1 -1  0  0;
#                         zeros(Ban, 1, 4) -1  1 -2 -2  0  2]];

A = [zeros(Ban,1,8) ones(Ban,1,2);
	 I M.*[-ones(Ban,4,1) zeros(Ban,4,1); zeros(Ban,4,1) -ones(Ban,4,1)]];

A_dom = [-1  1  1 -1 zeros(Ban, 1, 4)  0  0;
         zeros(Ban, 1, 4) -1  1  1 -1  0  0;
         zeros(Ban, 1, 4) -1  1 -2 -2  0  2];


A = [A; A_dom];
A = [A [zeros(Ban,1,11); I]] # only last 11 constraints are <=, the others are =

b = [1; zeros(Ban, 11)];

c, Q = load_param(experiment);

tol = 1e-8;
genLatex = false;
verbose = false;

sol = ipqp(A,b,c,Q, tol; maxit=15, verbose=verbose, genLatex=genLatex, slack_var=11:21);

print("\tSolution: ");
println([sol.x[1]-sol.x[2]+sol.x[5]-sol.x[6]; sol.x[3]-sol.x[4]+sol.x[7]-sol.x[8]]);
print("\tDisjoint flag: "); println(sol.x[9:10]);
println("")

nothing
