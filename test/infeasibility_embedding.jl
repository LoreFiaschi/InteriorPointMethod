include("../src/quadratic/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

num_var = 6;

Q = zeros(Ban, num_var, num_var);

c = [-1, -1, 0, 0, α, 0];
b = [-2, -1, -α];

A = [ 2 -1  1  0  -4  0;
     -1  2  0  1  -3  0;
	  0  0 -1 -1   0 -1];


A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);

tol=1e-8;
genLatex = false;
verbose = false;

if genLatex
	preamble();
end

sol = ipqp(A,b,c,Q, tol; maxit=15, verbose=verbose, genLatex=genLatex, slack_var=3:num_var);

if genLatex
	epilogue();
	println("");
	println("");
	println("");
	preamble();
	println("\t\\textbf{iter} & \$\\bm{r_1}\$ & \$\\bm{r_2}\$ & \$\\bm{r_3}\$ \\\\");
	println("\t\\hline");
	iter = 0
	for (r1, r2, r3) in eachrow(sol.r)
		global iter+=1;
		print("\t $(iter) & \$"); print_latex(r1); print("\$ & \$"); print_latex(r2); print("\$ & \$"); print_latex(r3); println("\$ \\\\"); 
		println("\t\\hline");
	end
	epilogue();
end
