function solve_standardlp(A,b,c,maxit=100,tol=1e-8,verbose=false,genLatex=false;adj=true,data_latex=[],slack_var=[])
    ### definition of gamma_f
	
	if genLatex && data_latex==[]
		throw(ArgumentError("to generate latex-wise output data_latex are needed."))
	end

    gamma_f = .01

	### definition of alpha_threshold
	
	alpha_threshold = (eltype(A)<:Real) ? 1e308 : α*α

    ### whether scaling is used
    
    scaling = 0

    m,n = size(A)
	
	var_to_show = setdiff(1:n, slack_var);

    ### compute initial value (x^0,lambda^0,s^0)
    
    x0,lambda0,s0 = starting_point(A,b,c)

    iter = 0
	
	r = Matrix(undef, 0, 3);

	#c_ref = deepcopy(c);
	c_inf = retrieve_infinitesimals(c, -1);
	c_fin = c-c_inf;

	if genLatex
		println("\t\\textbf{iter} & \$\\bm{\\mu}\$ & \\textbf{residual} & \$\\bm{x}\$ & \$\\bm{c^Tx}\$\\\\");
		println("\t\\hline");
		tmp_x = revProb(data_latex.r, get_x(data_latex.g, x0));
		print("\t$(iter) & \$"); print_latex(mean(x0.*s0)); print("\$ & \$"); print_latex(norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c])); print("\$ & \$"); print_latex(tmp_x[var_to_show]); print("\$ & \$"); print_latex(dot(data_latex.r.P.c, tmp_x)); println("\$ \\\\");
		println("\t\\hline");
    elseif verbose
        print(iter); print(" "); print(mean(x0.*s0)); print(" "); print(norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c])); print(" "); println("0., 0."); 
    end

    for iter=1:maxit
        ### solve 10.1
		
        f3 = fact3(A,x0,s0)

        # @show iter,ind_skip

        rb  = denoise(A*x0-b, tol)
        rc_fin = -c_fin
		rc_inf = denoise(A'*lambda0+s0-c_inf, tol)
        rxs = denoise(x0.*s0 , tol)
		#rxs = map(x->(x.p<-1) ? zero(Ban) : x,rxs)
		
		mu = mean(rxs)
		
		#=
		print("rb: "); println(rb);
		print("rc_inf: "); println(rc_inf);
		print("rxs: "); println(rxs);
		print("mu: "); println(mu)
		println("")
		=#

		#if iter >= 3
		#	lambda_aff_fin = zeros(m)
		#	x_aff_fin = s_aff_fin = zeros(n)
		#else
	        lambda_aff_fin,x_aff_fin,s_aff_fin = solve3(f3,A,x0,s0,zeros(m),rc_fin,zeros(n))
		#end
        lambda_aff_inf,x_aff_inf,s_aff_inf = solve3(f3,A,x0,s0,rb,rc_inf,rxs)

		#=
		println("")
		print("x_aff_fin: "); println(x_aff_fin)
		print("x_aff_inf: "); println(x_aff_inf)
		println("")
		=#	
		

		lambda_aff = lambda_aff_fin + lambda_aff_inf
		x_aff = x_aff_fin + x_aff_inf 
		s_aff = s_aff_fin + s_aff_inf

		#if iter>= 3
		#	lambda_cc = zeros(m)
		#	x_cc = s_cc = zeros(n)
		#else
		    ### calculate alpha_aff^pri, alpha_aff^dual, mu_aff

			#println("alpha_max x");
		    alpha_aff_pri  = alpha_max(x0,x_aff,1.0)
			
			#=
			println("");
			print("alpha_aff_pri: "); println(alpha_aff_pri)
			=#
		    alpha_aff_dual = alpha_max(s0,s_aff,1.0)

		    mu_aff = denoise(dot(x0+alpha_aff_pri*x_aff,s0+alpha_aff_dual*s_aff)/n , tol)
			#=
			print("mu_aff: "); println(mu_aff)
			println("")
			=#

		    ### centering parameter sigma

		    sigma = (mu_aff/mu)^3
		    
		    ### solve 10.7 

		    rb = zeros(m)
		    rc = zeros(n)
		    rxs = denoise(x_aff.*s_aff, tol).-sigma*mu 
			#rxs = map(x->(x.p<-1) ? zero(Ban) : x,rxs)

			#print("rxs :"); println(rxs)
			#println("")
			# print("rxs_aff: "); println(rxs);

		    lambda_cc,x_cc,s_cc = solve3(f3,A,x0,s0,rb,rc,rxs)
			#x_cc = map(x->(x.p<-1) ? zero(Ban) : x, x_cc)
			#lambda_cc = map(x->(x.p<-1) ? zero(Ban) : x, lambda_cc)
			#s_cc = map(x->(x.p<-1) ? zero(Ban) : x, s_cc)
			
			#=			
			print("x_cc: "); println(x_cc);
			print("s_cc: "); println(s_cc);
			=#
		#end
        ### compute search direction and step to boundary

        dx = x_aff+x_cc
        dlambda = lambda_aff+lambda_cc
        ds = s_aff+s_cc

        alpha_max_pri = alpha_max(x0,dx,Inf)		

        alpha_max_dual = alpha_max(s0,ds,Inf)

		#=
		print("alpha_aff_pri: "); println(alpha_aff_pri);
		print("alpha_aff_dual: "); println(alpha_aff_dual);
		print("alpha_max_pri: "); println(alpha_max_pri);
		print("alpha_max_dual: "); println(alpha_max_dual);
        =#
		
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
		
		
		#print("alpha_x: "); println(alpha_pri);
		#print("alpha_s: "); println(alpha_dual);	
		print("dx: "); println(dx);
		#print("ds: "); println(ds);
		println("");
		#println("");

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
		#print("x: "); println(x1[1:2])
        if genLatex
			tmp_x = revProb(data_latex.r, get_x(data_latex.g, x1));
			print("\t$(iter) & \$"); print_latex(mu); print("\$ & \$"); print_latex(norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c])); print("\$ & \$"); print_latex(tmp_x[var_to_show]); print("\$ & \$"); print_latex(dot(data_latex.r.P.c, tmp_x)); println("\$ \\\\");
			println("\t\\hline");
		elseif verbose
            # @printf("%3d %9.2e %9.2e %9.4g %9.4g\n", iter, mu, norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c]), alpha_pri, alpha_dual);
            print(iter); print(" "); print(mu); print(" "); print(norm([A'*lambda0 + s0 - c; A*x0 - b; x0.*s0])/norm([b;c])); print(" "); print(alpha_pri); print(" "); println(alpha_dual); 
        end

        ### termination

		#println(A*x1-b)
		r1 = denoise(norm(A*x1-b), tol)/(1+norm(b))
		r2 = denoise(norm(A'*lambda1+s1-c), tol)/(1+norm(c))
		cx = dot(c,x1)
		r3 = denoise(abs(cx-dot(b,lambda1)), tol)/(1+abs(cx))
		
		if genLatex			
			r = [r
				 r1 r2 r3];
	    end
        # @show r1
		
		
		#print("cx: "); println(cx);
		#print("bl: "); println(dot(b,lambda1));
		
		println("")
		print("r1: "); println(r1);
		print("r2: "); println(r2);
		print("r3: "); println(r3);
		println("");
		#=
		r1_d = denoise(r1,1e-15)
		if r1_d == 0 || r1_d.p < 0
			c = retrieve_infinitesimals(c_ref, -1)
			print("c: "); println(c)
		end
		println("")
		=#
		#end


        if (typeof(r1)<:Real) ? r1 < tol : all(z->abs(z) < tol, r1.num) 
		
            #r2 = norm(A'*lambda1+s1-c)/(1+norm(c))
            # @show r2
			
            if (typeof(r2)<:Real) ? r2 < tol : all(z->abs(z) < tol, r2.num) 

                #cx = dot(c,x1)
                #r3 = abs(cx-dot(b,lambda1))/(1+abs(cx))
                # @show r3

                if (typeof(r3)<:Real) ? r3 < tol : all(z->abs(z) < tol, r3.num) 
                    #println("");

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

    return x1,lambda1,s1,flag,iter,r
end
