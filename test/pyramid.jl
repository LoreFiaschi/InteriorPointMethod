include("../src/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

Q = [10 -2  4  0  0  0  0;
     -2 10  4  0  0  0  0;
	  4  4  4  0  0  0  0;
	  0  0  0  0  0  0  0;
	  0  0  0  0  0  0  0;
	  0  0  0  0  0  0  0;
	  0  0  0  0  0  0  0];


#c = [-16-η, -16, -16-η, 0, 0, 0, 0, 0]; # converges to (1,1,1)
c = [-16-η, -16-η, -16, 0, 0, 0, 0]; # converges to (1.5,1.5,0)

b = [1, 1, 1, 3];

A = [-1  1  1  1  0  0  0;  # the four faces of the pyramid
	  1  1 -1  0 -1  0  0;
	  1 -1  1  0  0  1  0;
	  1  1  1  0  0  0  1];

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     

tol=1e-8;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

sol = ipqp(A,b,c,Q, tol; maxit=25, verbose=verbose, genLatex=genLatex, slack_var=4:8);

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
