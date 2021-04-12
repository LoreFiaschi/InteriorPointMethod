include("../src/quadratic/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

function constraints3(W, H)
	# horizontal positioning
	A_h = [1  0 -1  0  0  0;
		   0  1 -1  0  0  0;
		   0  0  1  0  0  0];
	 
	b_h = [0; 0; W];


	# vertical positioning
	A_v = [0  0  0 -1  1  0;
		   0  0  0  1  0  0;
		   0  0  0  0  0  1];
	 
	b_v = [0; 10; H];


	# upper aspect ratio
	A_u = [2  0 -2 -1  0  0;
		   0  2 -2  1 -1  0;
		   0  0  2  0  0 -1];
	   
	b_u = [-H; 0; 2*W-H];


	# lower aspect ratio
	#=
	A_l = [0  0  0  1  0  0  0  0  0 -2  0  0;
		   0  0  0  0  1  0  0  0  0  0 -2  0;
		   0  0  0  0  0  1  0  0  0  0  0 -2];
	=#
	A_l = [-1  0  1  2  0  0;
		    0 -1  1 -2  2  0;
		    0  0 -1  0  0  2];

	b_l = [2*H; 0; 2*H-W];


	# minimal height
	A_m = [0  0  0  1  0  0;
		   0  0  0 -1  1  0;
		   0  0  0  0  0  1];
	   
	b_m = [H-1; -1; H-1];

	A = [A_h; A_v; A_u; A_l; A_m];
	
	b = [b_h; b_v; b_u; b_l; b_m];

	return A, b
end

function cost_f3(W, H)
	Q = -[0  0  0  1  0  0;
	      0  0  0 -1  1  0;
		  0  0  0  0 -1  1;
		  1 -1  0  0  0  0;
		  0  1 -1  0  0  0;
		  0  0  1  0  0  0]
	Q+= [ 4 -2 -4  0  0  0;
	     -2  2  0  0  0  0;
		 -4  0  8  0  0  0;
		  0  0  0 10 -4 -2;
		  0  0  0 -4  2  0;
		  0  0  0 -2  0  2].*η;
	
	c = -[-H, 0, 0, 0, 0, -W];
	c+=  [2*W, 0, -4*W, -4*H, 2*H, 0]*η
	
	return Q, c
end

function test(arr)
	for A in arr
		show(stdout, "text/plain", A);
		println("")
	end
end

function benchmark_3()
	H = 10
	W = 14
	
	A, b = constraints3(W,H)
	num_var = size(A,2)
	A = [A I]
	
	num_slack = size(A,1)

	Q, c = cost_f3(W,H)
	
	#=
	# Expected result: -276η
	x = [0, 0, 7, 4, 0, 0]
	print(0.5*x'*Q*x+c'*x)
	error()
	=#
	
	c = vcat(c, zeros(num_slack))
	Q = [Q zeros(size(Q,1), num_slack); zeros(num_slack, size(Q,1)+num_slack)]
	
	return A, b, c, Q, num_var
end

function benchmark_3_shift()
	H = 11
	W = 15
	
	A, b = constraints3(W,H)
	num_var = size(A,2)
	
	#shift
	A = [A; -I]
	b = vcat(b, -ones(num_var))
	
	A = [A I]
	num_slack = size(A,1)

	Q, c = cost_f3(W,H)
	
	#=
	# Expected result: 25-326η
	x = [1, 1, 8, 5, 1, 1]
	print(0.5*x'*Q*x+c'*x)
	error()
	=#
	
	c = vcat(c, zeros(num_slack))
	Q = [Q zeros(size(Q,1), num_slack); zeros(num_slack, size(Q,1)+num_slack)]
	
	return A, b, c, Q, num_var
end

#A, b, c, Q, num_var = benchmark_3()
A, b, c, Q, num_var = benchmark_3_shift()

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