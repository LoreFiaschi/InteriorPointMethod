include("../src/quadratic/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

using .BAN

# NOTICE!! Before launching assure that the tolerance considers the right powers for the stop criterion

#=
Q = [-1 -η 0 0 0 0
	 -η -1 0 0 0 0
	  ]
=#
  
Q = zeros(Ban, 6, 6);
Q[1,1] = -2;
Q[2,2] = -2;
#Q[1,2] = -η;
#Q[2,1] = -η;

c = [2-η, 2-η-η*η, 0, 0, 0, 0];
#c = [2, 2-η, 0, 0, 0, 0];
b = [1, 1, 3, -1];

A = [-1 1 1  0 0  0;    # y <=  x + 1
      1 1 0 -1 0  0;    # y >= -x + 1
      1 1 0  0 1  0;    # y <= -x + 3
     -1 1 0  0 0 -1];   # y >=  x - 1

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     

tol=1e-16;
genLatex = true;
verbose = false;

if genLatex
	preamble();
end

sol = ipqp(A,b,-c,-Q, tol; maxit=100, verbose=verbose, genLatex=genLatex, slack_var=3:6);

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
