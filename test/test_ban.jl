include("../src/iplp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

# NOTICE!! Before launching assure that the tolerance considers the right powers for the stop criterion

c = [-8-4η, -12-10η, 0, 0, 0, 0];
b = [120, 210, 270, 60];

A = [2 1 1 0 0  0;
     2 3 0 1 0  0;
     4 3 0 0 1  0;
     1 2 0 0 0 -1];

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     
     
lo = zeros(6);
#hi = [60;70;b[1:end-1];140];
hi = [Inf;Inf;Inf;Inf;Inf;Inf];

Problem = IplpProblem(c, A, b, lo, hi);

tol=1e-16;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

sol = iplp(Problem, tol; maxit=100, verbose=verbose, genLatex=genLatex);

if genLatex
	epilogue();
	println("");
	println("");
	println("");
	preamble();
	println("\t\$\\bm{r_1}\$ & \$\\bm{r_1}\$ & \$\\bm{r_1}\$ \\\\");
	println("\t\\hline");
	for (r1, r2, r3) in eachrow(sol.r)
		print("\t \$"); print_latex(r1); print("\$ & \$"); print_latex(r2); print("\$ & \$"); print_latex(r3); println("\$ \\\\"); 
		println("\t\\hline");
	end
	epilogue();
end