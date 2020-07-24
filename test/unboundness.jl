include("../src/iplp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using .BAN

# NOTICE!! Before launching assure that the tolerance considers the right powers for the stop criterion

c = [-1, -1, 0, 0];
b = [1, 1];

A = [-2  1 1  0;
      1 -1 0  1];

A = convert(Matrix{Ban}, A);     
#A = convert(SparseMatrixCSC{Ban}, A);     
     
lo = zeros(4);
hi = [α;α;2α+1;α+1];
#hi = [100;100;201;101];
#hi = [2;2;5;3];

Problem = IplpProblem(c, A, b, lo, hi);

tol=1e-8;
sol = iplp(Problem, tol; maxit=100);

println(sol.flag)
println(sol.x)