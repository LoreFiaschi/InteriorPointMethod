"""
Phase I to check the feasibility
"""

function phaseone(A, b, maxit=100, tol=1e-8)
    m,n = size(A)
    A = [A Matrix{Float64}(I,m,m)]
    c = [zeros(Float64, n);ones(Float64, m)]
    x1,lambda1,s1,flag,iter = solve_standardlp(A,b,c,maxit,tol,adj = false)
    # @show dot(c, x1)
    if abs(dot(c, denoise(x1, tol))) > 1e-8
		println(x1)
        return true
    end
    return false
end
