function solve_standardqp(A,b,c,Q, tol=1e-8, maxit=100; verbose=false, genLatex=false, adj=true, slack_var=[])
    ### definition of gamma_f

    gamma_f = .01

	### definition of α_threshold
	
	α_threshold = (eltype(A)<:Real) ? 1e308 : α*α

    m,n = size(A)
	
	var_to_show = setdiff(1:n, slack_var);
	
	flag = false;

	#################
    # initial value #
	#################
    
    x,λ,s = starting_point(A,b,c,Q)

    iter = 0
	
	r = Matrix(undef, 0, 3);

	if genLatex
		println("\t\\textbf{iter} & \$\\bm{\\μ}\$ & \\textbf{residual} & \$\\bm{x}\$ & \$\\bm{c^Tx}\$\\\\");
		println("\t\\hline");
		print("\t$(iter) & \$"); print_latex(mean(x.*s)); print("\$ & \$"); print_latex(norm([A'*λ + s - c; A*x - b; x.*s])/norm([b;c])); print("\$ & \$"); print_latex(x[var_to_show]); print("\$ & \$"); print_latex(0.5*(x'*Q*x)+dot(c, x)); println("\$ \\\\");
		println("\t\\hline");
    elseif verbose
        print(iter); print(" "); print(mean(x.*s)); print(" "); print(norm([A'*λ + s - c; A*x - b; x.*s])/norm([b;c])); print(" "); println("0., 0."); 
    end

    for iter=1:maxit
	
        ##############
		# solve 10.1 #
		##############
		
        f3 = fact3(A,Q,x,s)

        rb  = A*x-b 
        rc  = A'*λ+s-c-Q*x
        rxs = x.*s 

		μ = mean(rxs)

        λ_aff,x_aff,s_aff = solve3(f3,x,s,rb,rc,rxs)
		
		###########################
        # calculate α_aff, μ_aff #
		###########################

        α_aff_pri  = alpha_max(x,x_aff,1.0)
        α_aff_dual = alpha_max(s,s_aff,1.0)
		
		print("α_aff_pri: "); println(α_aff_pri); println("");
		print("α_aff_dual: "); println(α_aff_dual); println("");
		
        μ_aff = denoise(dot(x+α_aff_pri*x_aff,s+α_aff_dual*s_aff)/n , tol)
        σ = (μ_aff/μ)^3
		
		print("μ_aff: "); println(μ_aff); println("");
		print("σ: "); println(σ); println("");
        
		##############
        # solve 10.7 #
		##############

        rb = zeros(m)
        rc = zeros(n)
        rxs = x_aff.*s_aff.-σ*μ 

        λ_cc,x_cc,s_cc = solve3(f3,x,s,rb,rc,rxs)
		
		##############################
        # compute direction and step #
		##############################

        dx = x_aff+x_cc
        dλ = λ_aff+λ_cc
        ds = s_aff+s_cc
		
        α_pri = min(0.99*alpha_max(x,dx,Inf),1)
        α_dual = min(0.99*alpha_max(s,ds,Inf),1)
		
		print("α_pri: "); println(α_pri); println("");
		print("α_dual: "); println(α_dual); println("");

		#######################
        # check unboundedness #
		#######################

        if α_pri > α_threshold || α_dual > α_threshold
            @warn("This problem is unbounded")
            println(α_pri)
            println(α_dual)
            return x,λ,s,false,iter
        end

		###############################
        # compute x^k+1, λ^k+1, s^k+1 #
		###############################

        x = x+α_pri*dx
        λ = λ+α_dual*dλ
        s = s+α_dual*ds
		
		###############
        # termination #
		###############

		r1 = denoise(norm(A*x-b), tol)/(1+norm(b))
		r2 = denoise(norm(A'*λ+s-c-Q*x), tol)/(1+norm(c))
		r3 = denoise(dot(x,s), tol)/(1+abs(dot(c,x)+0.5*x'*Q*x))

        if genLatex
			print("\t$(iter) & \$"); print_latex(μ); print("\$ & \$"); print_latex(norm([A'*λ + s - c - Q*x; A*x - b; x.*s])/norm([b;c])); print("\$ & \$"); print_latex(x[var_to_show]); print("\$ & \$"); print_latex(0.5*(x'*Q*x)+dot(c, x)); println("\$ \\\\");
			println("\t\\hline");
			r = [r
				 r1 r2 r3];
		elseif verbose
            print(iter); print(" "); print(μ); print(" "); print(norm([A'*λ + s - c - Q*x; A*x - b; x.*s])/norm([b;c])); print(" "); print(α_pri); print(" "); println(α_dual); 
        end


        if (typeof(r1)<:Real) ? r1 < tol : all(z->abs(z) < tol, r1.num) 
		
            #r2 = norm(A'*λ+s-c)/(1+norm(c))			

            if (typeof(r2)<:Real) ? r2 < tol : all(z->abs(z) < tol, r2.num) 

                #cx = dot(c,x)
                #r3 = abs(cx-dot(b,λ))/(1+abs(cx))

                if (typeof(r3)<:Real) ? r3 < tol : all(z->abs(z) < tol, r3.num) 

					flag = true;
					break;
                end
            end
        end
    end

    return x,λ,s,flag,iter,r
end
