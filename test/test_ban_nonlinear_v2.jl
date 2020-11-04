include("../src/iplp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

# NOTICE!! Before launching assure that the tolerance considers the right powers for the stop criterion

#c = [2-η, 2-η-η*η, 0, 0, 0, 0];
#c = [-1, -1-η, 0, 0, 0, 0]; # converge to (1,2)
c = [-1-η, -1, 0, 0, 0, 0]; # converge to (2,1)
#c = [2, 2-η, 0, 0, 0, 0]; # converge to (0,1)
b = [1, 1, 3, -1];

A = [-1 1 1  0 0  0;    # y <=  x + 1
      1 1 0 -1 0  0;    # y >= -x + 1
      1 1 0  0 1  0;    # y <= -x + 3
     -1 1 0  0 0 -1];   # y >=  x - 1

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     

lo = zeros(6);
hi = [Inf;Inf;Inf;Inf;Inf;Inf];

Problem = IplpProblem(c, A, b, lo, hi);

tol=1e-8;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

sol = iplp(Problem, tol; maxit=100, verbose=verbose, genLatex=genLatex, slack_var=3:6);

if genLatex
	epilogue();
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
end
