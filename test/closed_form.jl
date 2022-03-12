include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s3_isbits.jl")
using .BAN
using LinearAlgebra

b = [120, 210, 270, 60];

A = [ 2  1;
	  2  3;
	  4  3;
	 -1 -2];

A = [A I]

A = convert(Matrix{Ban}, A);

n = size(A,2)

Q = zeros(Ban, n,n)
Q[1,1] = 1;
Q[2,2] = 1;

c = zeros(Ban, n);
c[1] = -40;
c[2] = -30;

M = A'*A+η*Q;
N = A'*b-η*c;

f = lu(M);
f = LU{eltype(M),typeof(M)}(denoise(f.factors, 1e-9), f.ipiv, f.info)
x = f\N;

println(x)
