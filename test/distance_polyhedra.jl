include("../src/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

function embed_2D(A_x, A_y, type="point";x=nothing,y=nothing)
	
	s_A_x = size(A_x)
	s_A_y = size(A_y)
	num_slack = s_A_x[1] + s_A_y[1];
	A = [A_x zeros(s_A_x[1], s_A_y[2]); zeros(s_A_y[1], s_A_x[2]) A_y];
	A = [A I];
	
	if type=="point"
		(x==nothing || y==nothing) && throw(ArgumentError("x and y must be specified in type point"))
		
		c = vcat(-2η*vcat(x,y), zeros(num_slack))
		Q = 2*[1+η  0  -1   0;
		        0  1+η  0  -1;
			   -1   0  1+η  0;
			    0  -1   0  1+η];
		Q = [Q zeros(4,num_slack);
		     zeros(num_slack, 4+num_slack)];
	elseif type=="barycenter"
		throw(ArgumentError("Not implemented yet"))
	end

	return A, Q, c
end

A_x = [-1  0;
       -1  1;
	    1  1;
		1 -1;
		0 -1];

A_y = [-3 -1;
	    1  3;
		1  1;
		3  1;
	   -1 -3;
	   -1 -1];

b = [-1, 3, 10, 2, -1.5, -25, 35, 17, 40, -20, -13];

x_ref = [2,5];
y_ref = [10.5, 4.5];

A, Q, c = embed_2D(A_x, A_y, "point", x=x_ref, y=y_ref)

tol=1e-8;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

sol = ipqp(A,b,c,Q, tol; maxit=20, verbose=verbose, genLatex=genLatex, slack_var=5:size(A,2));

if genLatex
	epilogue();
	#=
	println("");
	println("");
	println("");
	preamble();
	println("\t\\textbf{iter} & \$\\bm{r_1}\$ & \$\\bm{r_2}\$ & \$\\bm{r_3}\$ \\\\");
	println("\t\\hline");
	iter = 0;
	for (r1, r2, r3) in eachrow(sol.r)
		global iter+=1;
		print("\t $(iter) & \$"); print_latex(r1); print("\$ & \$"); print_latex(r2); print("\$ & \$"); print_latex(r3); println("\$ \\\\"); 
		println("\t\\hline");
	end
	epilogue();
	=#
end

nothing
