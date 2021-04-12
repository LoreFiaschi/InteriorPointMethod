include("../src/quadratic/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using LinearAlgebra
using .BAN

# useless
function box_constraints(num_cells, s)
	A = [I  I -ones(num_cells,1) zeros(num_cells, 2*num_cells+1);
	     zeros(num_cells, 2*num_cells+1) I  I -ones(num_cells,1)];
		 
	b = -ones(2*num_cells).*s
	
	return A, b
end

function horizontal_positioning(hp, s, num_cells, num_var)
	num_constraints = size(hp,1)
	A = zeros(num_constraints, num_var);
	for i=1:num_constraints
		A[i,hp[i,1]] = 1
		A[i,hp[i,1]+num_cells] = 1
		if hp[i,2]>num_cells
			A[i,num_cells*2+1] = -1
		else
			A[i,hp[i,2]] = -1
		end
	end
	b = -ones(num_constraints).*s
	
	return A, b
end

function vertical_positioning(vp, s, num_cells, num_var)
	num_constraints = size(vp,1)
	A = zeros(num_constraints, num_var);
	for i=1:num_constraints
		A[i,num_cells*2+1+vp[i,1]] = 1
		A[i,vp[i,1]+3*num_cells+1] = 1
		if vp[i,2]>num_cells
			A[i,end] = -1
		else
			A[i,num_cells*2+1+vp[i,2]] = -1
		end
	end
	b = -ones(num_constraints).*s
	
	return A, b
end

function aspect_ratio_upper(ar, num_cells, num_var)
	num_constraints = size(ar,1)
	A = zeros(num_constraints, num_var);
	for i=1:num_constraints
		A[i,convert(Int64, ar[i,1])+num_cells] = -ar[i,2]
		A[i,convert(Int64, ar[i,1])+num_cells*3+1] = 1
	end
	b = zeros(num_constraints)
	
	return A, b
end

function aspect_ratio_lower(ar, num_cells, num_var)
	num_constraints = size(ar,1)
	A = zeros(num_constraints, num_var);
	for i=1:num_constraints
		A[i,convert(Int64, ar[i,1])+num_cells] = ar[i,2]
		A[i,convert(Int64, ar[i,1])+num_cells*3+1] = -1
	end
	b = zeros(num_constraints)
	
	return A, b
end

function efficiency(num_cells, num_var)
	A = zeros(1, num_var)
	A[1,num_cells+1:2*num_cells] .= -1
	A[1,3*num_cells+2:end-1] .= -1
	b = -20
	
	return A, b
end

function minimal_height(num_cells, height)
	A = [zeros(num_cells, num_cells*3+1) -I zeros(num_cells, 1)]
	b = -ones(num_cells).*height
	
	return A, b
end

function perimeter(num_cells, num_var, num_slack)
	c = zeros(num_var+num_slack)
	c[num_cells*2+1] = 1
	c[num_cells*4+2] = 1
	
	return c
end

function area(num_cells, num_var, num_slack)
	Q = zeros(num_var+num_slack, num_var+num_slack)
	Q[num_cells*2+1, num_var] = 1
	Q[num_var, num_cells*2+1] = 1
	
	return Q
end

function area_coverage(num_cells, num_var, num_slack)
	Q = zeros(num_var+num_slack, num_var+num_slack)
	for i=1:num_cells
		Q[num_cells+i, 3*num_cells+1+i] = -1
		Q[3*num_cells+1+i, num_cells+i] = -1
	end
	
	return Q
end

function similarity(num_cells, num_var, num_slack)
	Q = zeros(num_var+num_slack, num_var+num_slack)
	
	# (N-1)*(w1^2+h1^2)
	Q[num_cells+1, num_cells+1] = 2*(num_cells-1)
	Q[num_cells*3+2, num_cells*3+2] = 2*(num_cells-1)
	
	for i=2:num_cells
		# wi^2+hi^2
		Q[num_cells+i, num_cells+i] = 2
		Q[num_cells*3+1+i, num_cells*3+1+i] = 2
		
		# -2w1*wi
		Q[num_cells+1,num_cells+i] = -2
		Q[num_cells+i,num_cells+1] = -2
		
		# -2h1*hi
		Q[num_cells*3+2,num_cells*3+1+i] = -2
		Q[num_cells*3+1+i,num_cells*3+2] = -2
	end
	
	return Q
end

#=
# START TEST

num_cells = 3
num_var = num_cells*4+2
s = 1

hp = [1 3;
	  2 3;
	  3 4];
	  
vp = [2 1;
	  1 4;
	  3 4];
	  
ar_u = [1 5;
        2 5;
        3 5]

ar_l = [1 0.2;
        2 0.2;
        3 0.2]
		
A, b = horizontal_positioning(hp, s, num_cells, num_var)

_A, _b = vertical_positioning(vp, s, num_cells, num_var)
A = [A; _A]
b = vcat(b,_b)

_A, _b = aspect_ratio_upper(ar_u, num_cells, num_var)
A = [A; _A]
b = vcat(b,_b)

_A, _b = aspect_ratio_lower(ar_l, num_cells, num_var)
A = [A; _A]
b = vcat(b,_b)

_A, _b = efficiency(num_cells, num_var)
A = [A; _A]
b = vcat(b,_b)

show(stdout, "text/plain", A);
println("")
show(stdout, "text/plain", b);
println("")

A = [A I]

num_slack = size(A,1)
c = perimeter(num_cells, num_slack)
show(stdout, "text/plain", c)
println("")

Q = area(num_cells, num_slack)
show(stdout, "text/plain", Q)
println("")

Q = similarity(num_cells, num_slack)
show(stdout, "text/plain", Q)
println("")

error()

# END TEST
=#

function benchmark_1()

	# number of cells to place
	num_cells = 5
	num_var = num_cells*4+2

	# spacing
	s = 1;

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
	ar_u = [1 5;
		   2 5;
		   3 5;
		   4 5;
		   5 5]

	ar_l = [1 0.2;
			2 0.2;
			3 0.2;
			4 0.2;
			5 0.2]

	# A: [x, w, W, y, h, H]

	A, b = horizontal_positioning(hp, s, num_cells, num_var)

	_A, _b = vertical_positioning(vp, s, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_upper(ar_u, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = aspect_ratio_lower(ar_l, num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	_A, _b = efficiency(num_cells, num_var)
	A = [A; _A]
	b = vcat(b,_b)

	A = [A I]

	num_slack = size(A,1)

	c = zeros(num_var+num_slack)
	#c = perimeter(num_cells, num_var, num_slack)*η
	Q = area(num_cells, num_var, num_slack)
	Q += area_coverage(num_cells, num_var, num_slack)*η
	#Q += similarity(num_cells, num_var, num_slack)*η#*η
	
	return A, b, c, Q, num_var
end


function benchmark_2()
	# number of cells to place
	num_cells = 5
	num_var = num_cells*4+2

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

	A, b = horizontal_positioning(hp, s, num_cells, num_var)

	_A, _b = vertical_positioning(vp, s, num_cells, num_var)
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

	A = [A I]

	num_slack = size(A,1)
	
	# box size constraints
	_A = zeros(2,size(A,2))
	_A[1,num_cells*2+1] = 1
	_A[2,num_var] = 1
	_b = [W, H]
	
	A = [A; _A]
	b = vcat(b,_b)


	c = zeros(num_var+num_slack)
	#c = perimeter(num_cells, num_var, num_slack)*η
	#Q = area(num_cells, num_var, num_slack)
	Q = area_coverage(num_cells, num_var, num_slack)
	Q += similarity(num_cells, num_var, num_slack)*η#*η
	
	return A, b, c, Q, num_var
end

#A, b, c, Q, num_var = benchmark_1()
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