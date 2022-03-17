include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s3_isbits.jl")
#include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN
using LinearAlgebra

function compute_increment(A1, b1, A2, b2, c, Q)

	M = A1'*A1 + η*Q + η*η*I
	K = A2'*A2 + M
	f = lu(M)
	N = f\(A1'*b1 - η*c)
	J = A2'*b2 - (A2'*A2)*N
	f = lu(K)
	return f\J
end
