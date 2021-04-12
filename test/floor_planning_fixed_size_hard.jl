include("../src/quadratic/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

function constraints()
	# horizontal positioning
	A_h = [0  0  0  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0;
		   1  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
		   0  1 -1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0;
		   0  0  1  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0;
		   0  0  0  1 -1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0];
	 
	b_h = [14; 0; 0; 0; 0];


	# vertical positioning
	A_v = [0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1;
		   0  0  0  0  0  0  0  0  0  0 -1  1  0  0  0  0  1  0  0  0;
		   0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  1  0  0  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0  0;];
	 
	b_v = [10; 10; 0; 0; 0];


	# upper aspect ratio
	A_u = [0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0  0  0  0;
		   0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0  0  0;
		   0  0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0  0;
		   0  0  0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0;
		   0  0  0  0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1;];
	   
	b_u = [0; 0; 0; 0; 0];


	# lower aspect ratio
	#=
	A_l = [0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0  0  0  0;
		   0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0  0  0;
		   0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0  0;
		   0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0;
		   0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2;];
	=#
	A_l = [0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;
		   0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0  0  0;
		   0  0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0  0;
		   0  0  0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0;
		   0  0  0  0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1;];
 
	b_l = [0; 0; 0; 0; 0];


	# minimal height
	A_m = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1;];
	   
	b_m = [-1; -1; -1; -1; -1];

	A = [A_h; A_v; A_u; A_l; A_m];
	
	b = [b_h; b_v; b_u; b_l; b_m];

	return A, b
end

function area_coverage(num_cells, num_var, num_slack)
	Q = zeros(num_var+num_slack, num_var+num_slack)
	for i=1:num_cells
		Q[num_cells+i, 3*num_cells+i] = -1
		Q[3*num_cells+i, num_cells+i] = -1
	end
	
	return Q
end

function similarity(num_cells, num_var, num_slack)
	Q = zeros(num_var+num_slack, num_var+num_slack)
	
	# (N-1)*(w1^2+h1^2)
	Q[num_cells+1, num_cells+1] = 2*(num_cells-1)
	Q[num_cells*3+1, num_cells*3+1] = 2*(num_cells-1)
	
	for i=2:num_cells
		# wi^2+hi^2
		Q[num_cells+i, num_cells+i] = 2
		Q[num_cells*3+i, num_cells*3+i] = 2
		
		# -2w1*wi
		Q[num_cells+1,num_cells+i] = -2
		Q[num_cells+i,num_cells+1] = -2
		
		# -2h1*hi
		Q[num_cells*3+1,num_cells*3+i] = -2
		Q[num_cells*3+i,num_cells*3+1] = -2
	end
	
	return Q
end

function benchmark_2()
	# number of cells to place
	num_cells = 5
	num_var = num_cells*4

	A, b = constraints()

	A = [A I]

	num_slack = size(A,1)

	c = zeros(num_var+num_slack)
	Q = area_coverage(num_cells, num_var, num_slack)
	
	# test to verify the correctness of the Q: must output 0.0
	#=
	x = ones(num_var+num_slack);
	for i=1:num_var+num_slack
		x[i] = i;
	end
	print(x'*Q*x+2*sum(x[num_cells+1:2*num_cells].*x[3*num_cells+1:4*num_cells]))
	println()
	error()
	=#
	_Q = similarity(num_cells, num_var, num_slack)
	Q += _Q*η
	#Q += similarity(num_cells, num_var, num_slack)*η#*η
	
	return A, b, c, Q, num_var
end

A, b, c, Q, num_var = benchmark_2()

tol=1e-8;
genLatex = false;
verbose = false;

if genLatex
	preamble();
end

sol = ipqp(A,b,c,Q, tol; maxit=50, verbose=verbose, genLatex=genLatex, slack_var=num_var+1:size(A,2));

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