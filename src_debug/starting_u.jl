"""
obtain the starting point for PD method
case of coupled variables
page 229
"""

function starting_point(A, b, c ,Q, bounded_variables, tol=1e-8)
    AA = denoise(A*A', tol)
    f = lu(AA)
	f = LU{eltype(AA),typeof(AA)}(denoise(f.factors, tol), f.ipiv, f.info)
	
	x = denoise(f\b, tol)
	x = denoise(A'*x, tol)
	
    λ = denoise(A*(c+Q*x), tol)
    λ = denoise(f\λ, tol)

    s = denoise(A[:,bounded_variables]'*λ, tol)
    s = denoise(c[bounded_variables]+Q[bounded_variables,:]*x-s, tol)
	
	print("true_x: ")
	println(x)
	println("")
	print("true_s: ")
	println(s)
	println("")
	

	# Version which cares about entries magnitude
	idx = findall(x->x<0, x)
	x[intersect(idx, bounded_variables)] .*= -0.5
	s[findall(x->x<0, s)] .*= -0.5

	xs = dot(x[bounded_variables],s)/2.0
	
    dx = xs/sum(s)
    ds = xs/sum(x[bounded_variables])

	mx = map(x->magnitude(x), x[bounded_variables])./magnitude(dx)

    x[bounded_variables] = x[bounded_variables]+(dx.*mx)
    s = s+(ds./mx)






	# Naive version
	#=
	# this may change entries magnitude
    dx = max(-1.5*minimum(x[bounded_variables]),0.0)
    ds = max(-1.5*minimum(s),0.0)
	x[bounded_variables] = x[bounded_variables].+dx
    s = s.+ds

    xs = dot(x[bounded_variables],s)/2.0
	
    dx = xs/sum(s)
    ds = xs/sum(x[bounded_variables])

	
	x[bounded_variables] = x[bounded_variables].+dx
    s = s.+ds
	=#


    return x,λ,s
end
