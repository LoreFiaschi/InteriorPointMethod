"""
obtain the starting point for PD method
(page 410 on Wright)
"""

function starting_point(A,b,c,Q)

    AA = A*A'
    f = cholesky(AA)

	x = f\b	
    x = A'*x
	
    λ = A*(c+Q*x)
    λ = f\λ

    s = A'*λ
    s = c+Q*x-s
	
    dx = max(-1.5*minimum(x),0.0)
    ds = max(-1.5*minimum(s),0.0)

    x = x.+dx
    s = s.+ds

    xs = dot(x,s)/2.0

    dx = xs/sum(s)
    ds = xs/sum(x)

    x = x.+dx
    s = s.+ds
	
	
    return x,λ,s
end
