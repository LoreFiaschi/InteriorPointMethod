include("../src/quadratic/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN
using LinearAlgebra

function linear_benchmark()
	Q = zeros(2,2)

	#c = [-8-14η, -12]; 	# converges to (30, 50)
	c = [-8-4η, -12-10η]; 	# converges to (0, 70)
	#c = [-8-14η, -12-10η]; 	# converges to (30, 50)
	#c = [-8, -12-10η]; 	# converges to (0, 70)

	return c, Q
end

function quadratic_benchmark()
	
	Q = [1 0; 0 1]*η # converges to (30, 50)
	#Q = [-1 0; 0 -1]*η # converges to (0, 70)
	
	c = [-8, -12]
	
	return c, Q
end

c, Q = linear_benchmark()
#c, Q = quadratic_benchmark()

b = [120, 210, 270, 60];

A = [ 2  1;
      2  3;
      4  3;
     -1 -2];

A = [A I]
c = vcat(c, zeros(size(A,1)))
Q = [Q zeros(size(Q,1), size(A,1)); zeros(size(A,1), size(Q,1)+size(A,1))]

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     

tol=1e-8;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);

if genLatex
	epilogue();
	println("");
	println("");
	println("");
	preamble();
	println("\t\$\\bm{r_1}\$ & \$\\bm{r_2}\$ & \$\\bm{r_3}\$ \\\\");
	println("\t\\hline");
	iter = 0;
	for (r1, r2, r3) in eachrow(sol.r)
		global iter+=1;
		print("\t $(iter) & \$"); print_latex(r1); print("\$ & \$"); print_latex(r2); print("\$ & \$"); print_latex(r3); println("\$ \\\\"); 
		println("\t\\hline");
	end
	epilogue();
end
