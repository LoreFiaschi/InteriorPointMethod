"""
obtain the starting point for PD method
(page 410 on Wright)
"""

function starting_point(A,b,c)
	#=
	println("");
	println("");
	print("A: "); println(A);
	println("");
	=#
    AA = A*A'

    f = cholesky(AA)
    # f = ldltfact(AA)
    # f = factorize(AA)

    # tilde
    x = f\b
    x = A'*x
	
    λ = A*c
    λ = f\λ

    s = A'*λ
    s = c-s
	#=
	print("x: "); println(x);
	println("");
	print("λ: "); println(λ);
	println("");
	print("s: "); println(s);
	println("");
	=#
    # hat
    dx = max(-1.5*minimum(x),0.0)
    ds = max(-1.5*minimum(s),0.0)

    x = x.+dx
    s = s.+ds

    # ^0
    xs = dot(x,s)/2.0

    dx = xs/sum(s)
    ds = xs/sum(x)
	#=
	print("dx: "); println(dx);
	println("");
	print("ds: "); println(ds);
	println("");
	println("");
	=#
    x = x.+dx
    s = s.+ds

    return x,λ,s
end
