# unreduced form

function fact3(A,Q,x,s)
    m,n = size(A)

    M = [zeros(m,m) A zeros(m,n);
	     A' -Q Matrix{Float64}(I,n,n);
	     zeros(n,m) diagm(s) diagm(x)]

    f = lu(M)

    return f
end

function solve3(f,x,s,rb,rc,rxs)
    m = length(rb)
    n = length(rc)

    b = Array([-rb; -rc; -rxs])

    b = f\b

    dλ = b[1:m]
    dx = b[1+m:m+n]
    ds = b[1+m+n:m+2*n]

    return dλ,dx,ds
end
