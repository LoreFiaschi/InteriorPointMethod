include("../src/iplp.jl")

c = [-8, -12, 0, 0, 0, 0];
b = [120, 210, 270, 60];

A = [2 1 1 0 0  0;
     2 3 0 1 0  0;
     4 3 0 0 1  0;
     1 2 0 0 0 -1];
     
A = convert(SparseMatrixCSC{Int64}, A);     
     
lo = zeros(6);
hi = [60;70;b]

Problem = IplpProblem(c, A, b, lo, hi);

tol=1e-8;
sol = iplp(Problem, tol; maxit=100);

println(sol.flag)
println(sol.x)
println(sol.iter)