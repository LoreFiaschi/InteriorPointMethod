function solve_standardlp(A,b,c,maxit=100,tol=1e-8,verbose=false;adj=true)
    ### definition of gamma_f

    gamma_f = .01
	
	### definition of alpha_threshold
	
	alpha_threshold = (eltype(A)<:Real) ? 1e308 : α

    ### whether scaling is used
    
    scaling = 1

    m,n = size(A)

    ### compute initial value (x^0,lambda^0,s^0)
    
    x0,lambda0,s0 = starting_point(A,b,c)

    # @show A*x0-b
    # @show s0

    iter = 0

    if verbose == true
        #@printf("%3d %9.2e %9.2e %9.4g %9.4g\n", iter, mean(x0.*s0), norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c]), 0., 0.)
        print(iter); print(" "); print(mean(x0.*s0)); print(" "); print(norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c])); print(" "); println("0., 0."); 
    end

    for iter=1:maxit
        ### solve 10.1

        f3 = fact3(A,x0,s0)

        # @show iter,ind_skip

        rb  = A*x0-b
        rc  = A'*lambda0+s0-c
        rxs = x0.*s0

        lambda_aff,x_aff,s_aff = solve3(f3,A,x0,s0,rb,rc,rxs)

        ### calculate alpha_aff^pri, alpha_aff^dual, mu_aff

        alpha_aff_pri  = alpha_max(x0,x_aff,1.0)
        alpha_aff_dual = alpha_max(s0,s_aff,1.0)

        mu = mean(rxs)
        mu_aff = dot(x0+alpha_aff_pri*x_aff,s0+alpha_aff_dual*s_aff)/n

        ### centering parameter sigma

        sigma = (mu_aff/mu)^3
        
        ### solve 10.7 

        rb = zeros(m)
        rc = zeros(n)
        rxs = x_aff.*s_aff.-sigma*mu

        lambda_cc,x_cc,s_cc = solve3(f3,A,x0,s0,rb,rc,rxs)

        ### compute search direction and step to boundary

        dx = x_aff+x_cc
        dlambda = lambda_aff+lambda_cc
        ds = s_aff+s_cc
		
		#print("dx: "); println(dx)

        alpha_max_pri = alpha_max(x0,dx,Inf)
        #=
        print("s0: "); println(s0);
        print("ds: "); println(ds);
        println("")
        =#
        alpha_max_dual = alpha_max(s0,ds,Inf)
        

        if scaling == 0
            alpha_pri = min(0.99*alpha_max_pri,1)
            alpha_dual = min(0.99*alpha_max_dual,1)
        else
            x1_pri = x0+alpha_max_pri*dx
            s1_dual = s0+alpha_max_dual*ds
            mu_p = dot(x1_pri,s1_dual)/n

            xind = argmin(x1_pri)
            # @test x1_pri[xind] == 0
            f_pri = (gamma_f*mu_p/s1_dual[xind]-x0[xind])/alpha_max_pri/dx[xind]
            sind = argmin(s1_dual)
            # @test s1_dual[sind] == 0
            f_dual = (gamma_f*mu_p/x1_pri[sind]-s0[sind])/alpha_max_dual/ds[sind]

            alpha_pri = max(1-gamma_f, f_pri)*alpha_max_pri
            alpha_dual = max(1-gamma_f, f_dual)*alpha_max_dual
        end

        # decide if the problem is unbounded

        if alpha_pri > alpha_threshold || alpha_dual > alpha_threshold
            @warn("This problem is unbounded")
            println(alpha_pri)
            println(alpha_dual)
            return x1,lambda1,s1,false,iter
        end

        ### compute x^k+1, lambda^k+1, s^k+1

        x1      = x0+alpha_pri*dx
        lambda1 = lambda0+alpha_dual*dlambda
        s1      = s0+alpha_dual*ds


        # @show dot(c,x1)
		#print("Obj: "); println(dot(c,x1))
		print("x: "); println(x1[1:2])
        #if verbose == true
            # @printf("%3d %9.2e %9.2e %9.4g %9.4g\n", iter, mu, norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c]), alpha_pri, alpha_dual);
            #print(iter); print(" "); print(mu); print(" "); print(norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c])); print(" "); print(alpha_pri); print(" "); println(alpha_dual); 
        #end

        ### termination

        r1 = norm(A*x1-b)/(1+norm(b))
        # @show r1
		print("r1: "); println(r1);

        if (typeof(r1)<:Real) ? r1 < tol : all(z->abs(z) < tol, r1.num[1:1+r1.p]) #r1 < tol #

            r2 = norm(A'*lambda1+s1-c)/(1+norm(c))
            # @show r2
			print("r2: "); println(r2);

            if (typeof(r2)<:Real) ? r2 < tol : all(z->abs(z) < tol, r2.num[1:1+r2.p]) #r2 < tol #

                cx = dot(c,x1)
                r3 = abs(cx-dot(b,lambda1))/(1+abs(cx))
                # @show r3
				print("r3: "); println(r3);

                if (typeof(r3)<:Real) ? r3 < tol : all(z->abs(z) < tol, r3.num[1:1+r3.p]) #r3 < tol #
                    println("");

					global flag = true
					global x1,lambda1,s1
					break
                end
            end
        end

        if iter == maxit
            global flag = false
            global x1,lambda1,s1
            break
        end

        x0      = x1
        lambda0 = lambda1
        s0      = s1

    end # end of for loop

    return x1,lambda1,s1,flag,iter
end
