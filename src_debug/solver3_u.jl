# unreduced form

function fact3(A,Q,x,s, bounded_variables)
    m,n = size(A)
	l_unbounded_variables = length(setdiff(1:n, bounded_variables))
	l_bounded_variables = length(bounded_variables)
	
#           dλ                              dz dx                                                                       ds
	M = [zeros(m,m)                           A                                                             zeros(m,l_bounded_variables);
	     A'                                  -Q                                                           [zeros(l_unbounded_variables, l_bounded_variables ); I]; 
	     zeros(l_bounded_variables,m) zeros(l_bounded_variables, l_unbounded_variables) diagm(s)                diagm(x[bounded_variables])]

    f = lu(M)
	f = LU{eltype(M),typeof(M)}(denoise(f.factors, tol), f.ipiv, f.info)
	
    return f
end

function solve3(f,rb,rc,rxs)
    m = length(rb)
    n = length(rc)

    b = f\[rb; rc; rxs]

    dλ = b[1:m]
    dx = b[1+m:m+n]
    ds = b[1+m+n:end]

    return dλ,dx,ds
end
