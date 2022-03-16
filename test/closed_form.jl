include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s3_isbits.jl")
#include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN
using LinearAlgebra

function load_params(experiment)

	if experiment == 0 # kite inner region, optimal point (40, 30)

		b = [120, 210, 270, -60];

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

		M = A'*A+η*Q

		return A, b, c, Q, M
	
	elseif experiment == 1 # kite outer region, optimal point (50, 50)

		b = [120, 210, 270, -60];

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
		c[1] = -50;
		c[2] = -50;


		M = A'*A+η*Q

		return A, b, c, Q, M

	elseif experiment == 2 # pg 452 example 16.2, optimal point (2, -1, 1)

		A = [1 0 1;
			 0 1 1];

		b = [3; 0];

		Q = [6 2 1;
			 2 5 2;
			 1 2 4];

		c = [-8; -3; -3];

		M = A'*A+η*Q

		return A, b, c, Q, M

	elseif experiment == 3

		throw(ArgumentError("Experiment $(experiment) not yet implemented (degenerate one)"))

	else
		throw(ArgumentError("Experiment $(experiment) does not exists"))
	end

end

experiment = 2



A, b, c, Q, M = load_params(experiment)

N = A'*b-η*c;

f = lu(M);
f = LU{eltype(M),typeof(M)}(denoise(f.factors, 1e-9), f.ipiv, f.info)
#f = cholesky(M)
x = f\N;

println("Solution:")
display(x)
println("")
println("Feasibility:")
display(A*x-b)
println("")
println("Value function:")
display(0.5*x'*Q*x+c'x)
println("")
