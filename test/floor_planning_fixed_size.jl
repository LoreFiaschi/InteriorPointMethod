include("../src/ipqp.jl")
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
		   0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0  0];
	 
	b_v = [10; 10; 0; 0; 0];


	# upper aspect ratio
	A_u = [0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0  0  0  0;
		   0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0  0  0;
		   0  0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0  0;
		   0  0  0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1  0;
		   0  0  0  0  0  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  1];
	   
	b_u = [0; 0; 0; 0; 0];


	# lower aspect ratio
	#=
	A_l = [0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0  0  0  0;
		   0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0  0  0;
		   0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0  0;
		   0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2  0;
		   0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -2];
	=#
	A_l = [0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;
		   0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0  0  0;
		   0  0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0  0;
		   0  0  0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1  0;
		   0  0  0  0  0  0  0  0  0 .5  0  0  0  0  0  0  0  0  0 -1]; 
	b_l = [0; 0; 0; 0; 0];


	# minimal height
	A_m = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;
		   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1];
	   
	b_m = [-1; -1; -1; -1; -1];

	A = [A_h; A_v; A_u; A_l; A_m];
	
	b = [b_h; b_v; b_u; b_l; b_m];

	return A, b
end

function constraints3()
	# horizontal positioning
	A_h = [0  0  1  0  0  1  0  0  0  0  0  0;
		   1  0 -1  1  0  0  0  0  0  0  0  0;
		   0  1 -1  0  1  0  0  0  0  0  0  0];
	 
	b_h = [14; 0; 0];


	# vertical positioning
	A_v = [0  0  0  0  0  0  1  0  0  1  0  0;
		   0  0  0  0  0  0  0  0  1  0  0  1;
		   0  0  0  0  0  0 -1  1  0  0  1  0];
	 
	b_v = [10; 10; 0];


	# upper aspect ratio
	A_u = [0  0  0 -2  0  0  0  0  0  1  0  0;
		   0  0  0  0 -2  0  0  0  0  0  1  0;
		   0  0  0  0  0 -2  0  0  0  0  0  1];
	   
	b_u = [0; 0; 0];


	# lower aspect ratio
	#=
	A_l = [0  0  0  1  0  0  0  0  0 -2  0  0;
		   0  0  0  0  1  0  0  0  0  0 -2  0;
		   0  0  0  0  0  1  0  0  0  0  0 -2];
	=#
	A_l = [0  0  0 .5  0  0  0  0  0 -1  0  0;
		   0  0  0  0 .5  0  0  0  0  0 -1  0;
		   0  0  0  0  0 .5  0  0  0  0  0 -1];

	b_l = [0; 0; 0];


	# minimal height
	A_m = [0  0  0  0  0  0  0  0  0 -1  0  0;
		   0  0  0  0  0  0  0  0  0  0 -1  0;
		   0  0  0  0  0  0  0  0  0  0  0 -1];
	   
	b_m = [-1; -1; -1];

	A = [A_h; A_v; A_u; A_l; A_m];
	
	b = [b_h; b_v; b_u; b_l; b_m];

	return A, b
end

function test(arr)
	for A in arr
		show(stdout, "text/plain", A);
		println("")
	end
end

function horizontal_positioning(hp, s, num_cells, num_var, W, bias=0)
	num_constraints = size(hp,1)
	A = zeros(num_constraints, num_var);
	b = (s==0) ? zeros(num_constraints) : -ones(num_constraints).*s
	W += bias
	
	for i=1:num_constraints
		A[i,hp[i,1]] = 1
		A[i,hp[i,1]+num_cells] = 1
		if hp[i,2]>num_cells
			b[i] += W
		else
			A[i,hp[i,2]] = -1
		end
	end
	
	if bias!=0
		A = [A; -I zeros(num_cells, num_cells*3)]
		b = vcat(b, -ones(num_cells).*bias)
	end
	
	return A, b
end

function vertical_positioning(vp, s, num_cells, num_var, H, bias=0)
	num_constraints = size(vp,1)
	A = zeros(num_constraints, num_var);
	b = (s==0) ? zeros(num_constraints) : -ones(num_constraints).*s
	H += bias
	
	for i=1:num_constraints
		A[i,num_cells*2+vp[i,1]] = 1
		A[i,vp[i,1]+3*num_cells] = 1
		if vp[i,2]>num_cells
			b[i] += H
		else
			A[i,num_cells*2+vp[i,2]] = -1
		end
	end
	
	if bias!=0
		A = [A; zeros(num_cells, num_cells*2) -I zeros(num_cells, num_cells)]
		b = vcat(b, -ones(num_cells).*bias)
	end
	
	return A, b
end

function aspect_ratio_upper(ar, num_cells, num_var)
	num_constraints = size(ar,1)
	A = zeros(num_constraints, num_var);
	for i=1:num_constraints
		A[i,convert(Int64, ar[i,1])+num_cells] = -ar[i,2]
		A[i,convert(Int64, ar[i,1])+num_cells*3] = 1
	end
	b = zeros(num_constraints)
	
	return A, b
end

function aspect_ratio_lower(ar, num_cells, num_var)
	num_constraints = size(ar,1)
	A = zeros(num_constraints, num_var);
	for i=1:num_constraints
		A[i,convert(Int64, ar[i,1])+num_cells] = ar[i,2]
		A[i,convert(Int64, ar[i,1])+num_cells*3] = -1
	end
	b = zeros(num_constraints)
	
	return A, b
end

function minimal_height(num_cells, height)
	A = [zeros(num_cells, num_cells*3) -I]
	b = -ones(num_cells).*height
	
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

	# spacing
	s = 0
	
	# minimal height
	min_h = 1
	
	# box size
	H = 10
	W = 14

	# horizontal positioning (to be on the left of)
	hp = [5 6;
		  1 3;
		  2 3;
		  3 5;
		  4 5];

	# vertical positioning (to be under to)
	vp = [4 6;
		  5 6;
		  2 1;
		  1 4;
		  3 4];
		 
	# aspect ratio
	ar_u = [1 2;
		   2 2;
		   3 2;
		   4 2;
		   5 2]

	ar_l = [1 0.5;
			2 0.5;
			3 0.5;
			4 0.5;
			5 0.5]

	# A: [x, w, W, y, h, H]

	A, b = horizontal_positioning(hp, s, num_cells, num_var, W)

	_A, _b = vertical_positioning(vp, s, num_cells, num_var, H)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_upper(ar_u, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_lower(ar_l, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = minimal_height(num_cells, min_h)
	A = [A; _A]
	b = vcat(b,_b)
	
	# test with hard-coded version: must output true
	#=
	_A, _ = constraints()
	print(A-_A == zeros(size(A)))
	println("")
	error()
	=#

	A = [A I]

	num_slack = size(A,1)

	c = zeros(num_var+num_slack)
	Q = area_coverage(num_cells, num_var, num_slack)
	_Q = similarity(num_cells, num_var, num_slack)
	Q += _Q*η
	#Q += similarity(num_cells, num_var, num_slack)*η#*η
	
	return A, b, c, Q, num_var
end

function benchmark_3()
	# number of cells to place
	num_cells = 3
	num_var = num_cells*4

	# spacing
	s = 0
	
	# minimal height
	min_h = 1
	
	# box size
	H = 10
	W = 14

	# horizontal positioning (to be on the left of)
	hp = [3 4;
		  1 3;
		  2 3];

	# vertical positioning (to be under to)
	vp = [1 4;
		  3 4;
		  2 1];
		 
	# aspect ratio
	ar_u = [1 2;
		    2 2;
		    3 2]

	ar_l = [1 0.5;
			2 0.5;
			3 0.5]

	# A: [x, w, W, y, h, H]

	A, b = horizontal_positioning(hp, s, num_cells, num_var, W)

	_A, _b = vertical_positioning(vp, s, num_cells, num_var, H)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_upper(ar_u, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_lower(ar_l, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = minimal_height(num_cells, min_h)
	A = [A; _A]
	b = vcat(b,_b)

	# test with hard-coded version: must output true
	#=
	_A, _ = constraints3()
	show(stdout, "text/plain", A)
	println("")
	show(stdout, "text/plain", _A)
	println("")
	print(A-_A == zeros(size(A)))
	println("")
	error()
	=#
	
	# test to verify feasibility of optimum and its optimality
	#=
	x0 = [0; 0; 3; 3; 3; 11; 5; 0; 0; 5; 5; 10];
	# must be non-negative
	show(stdout, "text/plain", b-A*x0)
	println("")
	# any(x->x<0, b-A*x0) ? println(false) : println(true)
	# must output -280
	c = zeros(Ban, num_var)
	Q = area_coverage(num_cells, num_var, 0)
	show(stdout, "text/plain", x0'*Q*x0+c'*x0)
	println("")
	error()
	=#
	
	A = [A I]
	

	num_slack = size(A,1)

	c = zeros(Ban, num_var+num_slack)
	Q = area_coverage(num_cells, num_var, num_slack)
	#_Q = similarity(num_cells, num_var, num_slack)
	#Q += _Q.*η
	#Q = convert(Array{Ban, 2}, Q)
	Q += similarity(num_cells, num_var, num_slack)*η#*η
	
	return A, b, c, Q, num_var
end

# starting from (1,1)
function benchmark_4()
	# initial gap
	bias = 1

	# number of cells to place
	num_cells = 3
	num_var = num_cells*4

	# spacing
	s = 0
	
	# minimal height
	min_h = 1
	
	# box size
	H = 10
	W = 14

	# horizontal positioning (to be on the left of)
	hp = [3 4;
		  1 3;
		  2 3];

	# vertical positioning (to be under to)
	vp = [1 4;
		  3 4;
		  2 1];
		 
	# aspect ratio
	ar_u = [1 2;
		    2 2;
		    3 2]

	ar_l = [1 0.5;
			2 0.5;
			3 0.5]

	# A: [x, w, W, y, h, H]

	A, b = horizontal_positioning(hp, s, num_cells, num_var, W, bias)

	_A, _b = vertical_positioning(vp, s, num_cells, num_var, H, bias)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_upper(ar_u, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_lower(ar_l, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = minimal_height(num_cells, min_h)
	A = [A; _A]
	b = vcat(b,_b)

	# test with hard-coded version: must output true
	#=
	_A, _ = constraints3()
	show(stdout, "text/plain", A)
	println("")
	show(stdout, "text/plain", _A)
	println("")
	print(A-_A == zeros(size(A)))
	println("")
	error()
	=#
	
	A = [A I]
	

	num_slack = size(A,1)

	c = zeros(Ban, num_var+num_slack)
	Q = area_coverage(num_cells, num_var, num_slack)
	#_Q = similarity(num_cells, num_var, num_slack)
	#Q += _Q.*η
	#Q = convert(Array{Ban, 2}, Q)
	Q += similarity(num_cells, num_var, num_slack)*η#*η
	
	return A, b, c, Q, num_var
end

#A, b, c, Q, num_var = benchmark_2()
A, b, c, Q, num_var = benchmark_3()
#A, b, c, Q, num_var = benchmark_4()

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
