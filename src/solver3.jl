# unreduced form

function fact3(A,x,s)
    m,n = size(A)

    M = [zeros(m,m) A zeros(m,n);
	     A' zeros(n,n) Matrix{Float64}(I,n,n);
	     zeros(n,m) diagm(s[:,1]) diagm(x[:,1])]

#=
	println(s);
	println("");
	println(diagm(0=>s[:,1]));
	println("");
	println(x);
	println("");
	println(diagm(0=>x[:,1]))
	println("");
	error();
=#
	#println(M);
	#println("");
    f = lu(M)
	#=
	println("")
	println("U")
	println(f.U)
	println("")
	println("L")
	println(f.L)
	println("")	
	println("")
	error()
	=#
    return f
end

function solve3(f,A,x,s,rb,rc,rxs)
    m = length(rb)
    n = length(rc)

    b = Array([-rb; -rc; -rxs])

    b = f\b

    dlam = b[1:m]
    dx = b[1+m:m+n]
    ds = b[1+m+n:m+2*n]

    return dlam,dx,ds
end
