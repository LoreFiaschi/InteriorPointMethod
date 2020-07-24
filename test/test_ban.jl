include("../src/iplp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

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
hi = [60;70;b[1:end-1];Inf];

Problem = IplpProblem(c, A, b, lo, hi);

tol=1e-8;
#tol = Ban(0, ones(3).*1e-8, false);
sol = iplp(Problem, tol; maxit=100);

println(sol.flag)
println(sol.x)