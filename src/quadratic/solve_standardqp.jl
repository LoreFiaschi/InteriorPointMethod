function solve_standardqp(A,b,c,Q, tol=1e-8, maxit=100; verbose=false, genLatex=false, adj=true, slack_var=[])
    ### definition of gamma_f

    gamma_f = .01

	### definition of α_threshold
	
	α_threshold = (eltype(A)<:Real) ? 1e308 : α*α

    m,n = size(A)
	trash_deg = minimum(map(x->min_degree(x), Q))-1 # to improve considering also b and A
	trash_deg = min(trash_deg, minimum(map(x->min_degree(x), c))-1) # to improve considering also b and A
	
	var_to_show = setdiff(1:n, slack_var);
	
	flag = false;

	#################
    # initial value #
	#################
    
    x,λ,s = starting_point(A,b,c,Q)
	
	print("x0: "); println(x)
	println("")
	#print("s0: "); println(s)
	#println("")
	#print("λ0: "); println(λ)
	#println("")

	# Cleaning
	
	x = denoise(x, tol)
	x[findall(x->x<0, x)] .= 0
	x -= retrieve_infinitesimals(x, trash_deg)
	
	s = denoise(s, tol)
	s[findall(x->x<0, s)] .= 0
	s -= retrieve_infinitesimals(s, trash_deg)
		
	λ = denoise(λ, tol)
	λ -= retrieve_infinitesimals(λ, trash_deg)

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

        #rb  = denoise(A*x-b, tol)
        #rc  = denoise(A'*λ+s-c-Q*x, tol)

		print("x: "); println(x)
		println("")
		#print("s: "); println(s)
		#println("")
		#print("λ: "); println(λ)
		#println("")

		
		rb  = denoise(b-A*x, tol)
        rc  = denoise(c+Q*x-A'*λ-s, tol)
        rxs = denoise(-x.*s, tol)

		#rb  = b-A*x
        #rc  = c+Q*x-A'*λ-s
        #rxs = -x.*s
		
		print("n_rb: "); println(norm(rb))
		print("n_rc: "); println(norm(rc))
		print("n_rxs: "); println(norm(rxs))
		println("")

		μ = mean(-rxs)
		
		print("μ: "); println(μ);
		println("")

        λ_aff,x_aff,s_aff = solve3(f3,rb,rc,rxs)
		
		
		print("n_x_aff_true: "); println(norm(x_aff))
		println("")
		#print("n_s_aff_true: "); println(norm(s_aff))
		#println("")
		#print("n_λ_aff_true: "); println(norm(λ_aff))
		#println("")
		
		x_aff = denoise(x_aff, tol) 
		s_aff = denoise(s_aff, tol)
		λ_aff = denoise(λ_aff, tol)		
		x_aff[findall(x->x.p<trash_deg, x_aff)] .= 0
		s_aff[findall(x->x.p<trash_deg, s_aff)] .= 0
		
		print("x_aff: "); println(x_aff)
		println("")
		#print("s_aff: "); println(s_aff)
		#println("")
		#print("λ_aff: "); println(λ_aff)
		#println("")

		###########################
        # calculate α_aff, μ_aff #
		###########################

        α_aff_pri  = alpha_max(x,x_aff,1.0)
        α_aff_dual = alpha_max(s,s_aff,1.0)
		
		print("α_aff_pri: "); println(α_aff_pri)
		print("α_aff_dual: "); println(α_aff_dual)
		println("")
				
        μ_aff = dot(denoise(x+α_aff_pri*x_aff, tol), denoise(s+α_aff_dual*s_aff, tol))/n
        (μ==0) ? σ = 0 : σ = (μ_aff/μ)^3 # σ = (μ_aff/μ)^3 # 
		
		print("μ_aff: "); println(μ_aff)
		print("σ: "); println(σ)
		println("")
				
		##############
        # solve 10.7 #
		##############

        rb = zeros(m)
        rc = zeros(n)
        rxs = σ*μ.-x_aff.*s_aff
		
		print("n_rxs: "); println(norm(rxs))
		println("")

        λ_cc,x_cc,s_cc = solve3(f3,rb,rc,rxs)
		
		print("n_x_cc_true: "); println(norm(x_cc))
		println("")
		#print("n_s_cc_true: "); println(norm(s_cc))
		#println("")
		#print("n_λ_cc_true: "); println(norm(λ_cc))
		#println("")

		x_cc = denoise(x_cc, tol) 
		s_cc = denoise(s_cc, tol)
		λ_cc = denoise(λ_cc, tol)		
		x_cc[findall(x->x.p<trash_deg, x_cc)] .= 0
		s_cc[findall(x->x.p<trash_deg, s_cc)] .= 0
		
		print("x_cc: "); println(x_cc)
		println("")
		#print("s_cc: "); println(s_cc)
		#println("")
		#print("λ_cc: "); println(λ_cc)
		#println("")
		
		##############################
        # compute direction and step #
		##############################

        dx = x_aff+x_cc
        dλ = λ_aff+λ_cc
        ds = s_aff+s_cc
		
		print("dx: "); println(dx)
		println("")
		#print("ds: "); println(ds)
		#println("")
		#print("dλ: "); println(dλ)
		#println("")
		
        α_pri = min(0.99*alpha_max(x,dx,Inf),1)
        α_dual = min(0.99*alpha_max(s,ds,Inf),1)
		
		print("α_pri: "); println(α_pri)
		print("α_dual: "); println(α_dual)
		println("")
		
		###############################
        # compute x^k+1, λ^k+1, s^k+1 #
		###############################

        x = x+α_pri*dx
        λ = λ+α_dual*dλ
        s = s+α_dual*ds
		
		x = denoise(x, tol)
		x[findall(x->x<0, x)] .= 0
		x -= retrieve_infinitesimals(x, trash_deg)
		
		s = denoise(s, tol)
		s[findall(x->x<0, s)] .= 0
		s -= retrieve_infinitesimals(s, trash_deg)
			
		λ = denoise(λ, tol)
		λ -= retrieve_infinitesimals(λ, trash_deg)
		
		###############
        # termination #
		###############

		r1 = denoise(norm(A*x-b), tol)/(1+norm(b))
		r2 = denoise(norm(A'*λ+s-c-Q*x), tol)/(1+norm(c))
		r3 = denoise(dot(x,s)/n, tol)/(1+abs(dot(c,x)+0.5*x'*Q*x))
		#r3 = denoise(abs(dot(c,x)-dot(b,λ)+x'*Q*x), tol)/(1+abs(dot(c,x)+0.5*x'*Q*x))
		
		r1 -= retrieve_infinitesimals(r1, trash_deg)
		r2 -= retrieve_infinitesimals(r2, trash_deg)
		r3 -= retrieve_infinitesimals(r3, trash_deg)

        if genLatex
			print("\t$(iter) & \$"); print_latex(μ); print("\$ & \$"); print_latex(norm([A'*λ + s - c - Q*x; A*x - b; x.*s])/norm([b;c])); print("\$ & \$"); print_latex(x[var_to_show]); print("\$ & \$"); print_latex(0.5*(x'*Q*x)+dot(c, x)); println("\$ \\\\");
			println("\t\\hline");
			r = [r
				 r1 r2 r3];
		elseif verbose
            print(iter); print(" "); print(μ); print(" "); print(norm([A'*λ + s - c - Q*x; A*x - b; x.*s])/norm([b;c])); print(" "); print(α_pri); print(" "); println(α_dual); 
        end

		#println("")
		#print("r1: "); println(r1)
		#print("r2: "); println(r2)
		#print("r3: "); println(r3)
		#println("");


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
