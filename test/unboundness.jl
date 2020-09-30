include("../src/iplp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

# NOTICE!! Before launching assure that the tolerance considers the right powers for the stop criterion

c = [-1, -1, 0, 0, 0];
#b = [1, 1, 1000];
b = [1, 1, α];
#b = [η, η, 1];

A = [-2  1 1  0  0;
      1 -1 0  1  0;
	  1  1 0  0  1];


	  
lo = zeros(5);
#hi = [α/2;α/2;α+1;α/2+1;α];
hi = [Inf;Inf;Inf;Inf;Inf];
#=
c = [-1, -1, 0, 0];
b = [1, 1];

A = [-2  1 1  0;
      1 -1 0  1];
	  
lo = zeros(4);
hi = [α;α;2α+1;α+1];
=#

#hi = [100;100;201;101];
#hi = [2;2;5;3];

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     

Problem = IplpProblem(c, A, b, lo, hi);

tol=1e-16;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

tol=1e-16;
sol = iplp(Problem, tol; maxit=100, verbose=verbose, genLatex=genLatex, slack_var=3:5);

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
