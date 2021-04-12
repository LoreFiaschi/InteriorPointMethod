function solve_standardqp(A,b,c,Q, tol=1e-8, maxit=100; verbose=false, genLatex=false, adj=true, slack_var=[])

	###########################
	# compute round tresholds #
	###########################

    m,n = size(A)
	min_deg_Q = minimum(map(x->min_degree(x), Q))
	min_deg_A = minimum(map(x->min_degree(x), A))
	min_deg_c = minimum(map(x->min_degree(x), c))
	min_deg_b = minimum(map(x->min_degree(x), b))
	
	trash_deg_r1 = min(min_deg_A, min_deg_b)-1;
	trash_deg_r2 = min(min_deg_Q, min_deg_A, min_deg_c)-1;
	trash_deg = min(trash_deg_r2, min_deg_b);
	
	#####################
	# garbage variables #
	#####################
	
	var_to_show = setdiff(1:n, slack_var);
	
	flag = false;

	#################
    # initial value #
	#################
    
    x,λ,s = starting_point(A,b,c,Q)

	############
	# Cleaning #
	############
	
	x = denoise(x, tol)
	x[findall(x->x<0, x)] .= 0
	x -= retrieve_infinitesimals(x, trash_deg)
	
	s = denoise(s, tol)
	s[findall(x->x<0, s)] .= 0
	s -= retrieve_infinitesimals(s, trash_deg)
		
	λ = denoise(λ, tol)
	λ -= retrieve_infinitesimals(λ, trash_deg)
	
	#################
	# Aux variables #
	#################

    iter = 0
	show = true
	show_more = false
	r = Matrix(undef, 0, 3); # just for genLatex purposes
	
	rb_den = norm(b);
	rb_den += magnitude(rb_den);
	rc_den = norm(c);
	rc_den += magnitude(rc_den);

	if genLatex
		println("\t\\textbf{iter} & \$\\bm{\\mu}\$ & \$\\bm{x}\$ & \$\\bm{c^Tx}\$\\\\");
		println("\t\\hline");
		print("\t$(iter) & \$"); print_latex(mean(x.*s)); print("\$ & \$"); print_latex(x[var_to_show]); print("\$ & \$"); print_latex(0.5*(x'*Q*x)+dot(c, x)); println("\$ \\\\");
		println("\t\\hline");
    elseif verbose
        print(iter); print(" "); print(mean(x.*s)); print(" "); print(norm([A'*λ + s - c; A*x - b; x.*s])/norm([b;c])); print(" "); println("0., 0."); 
    end

    for iter=1:maxit
	
        ##############
		# solve 10.1 #
		##############

		rb  = b-A*x
        rc  = c+Q*x-A'*λ-s
        rxs = -x.*s
		
		rb = denoise(rb, tol)
		rc = denoise(rc, tol)
		rxs = denoise(rxs, tol)
		#rb[findall(x->x.p<=trash_deg, rb)] .= 0
		#rc[findall(x->x.p<=trash_deg, rc)] .= 0
		#rxs[findall(x->x.p<=trash_deg, rxs)] .= 0
		rb -= retrieve_infinitesimals(rb, trash_deg_r1)
		rc -= retrieve_infinitesimals(rc, trash_deg_r2)
		rxs -= retrieve_infinitesimals(rxs, trash_deg)

		
		print("rb: "); println(rb)
		println("")
		print("rc: "); println(rc)
		println("")
		print("rxs: "); println(rxs)
		println("")
		
		if show_more
			print("rb: "); println(norm(rb))
			println("")
			print("rc: "); println(norm(rc))
			println("")
		end
		
		
	
		f3 = fact3(A,Q,x,s)
        λ_aff,x_aff,s_aff = solve3(f3,rb,rc,rxs)
		#=
		print("x_aff: "); println(x_aff)
		println("")
		print("s_aff: "); println(s_aff)
		println("")
		print("λ_aff: "); println(λ_aff)
		println("")
		=#
		x_aff = denoise(x_aff, tol) 
		s_aff = denoise(s_aff, tol)
		λ_aff = denoise(λ_aff, tol)		
		x_aff[findall(x->x.p<=trash_deg, x_aff)] .= 0
		s_aff[findall(x->x.p<=trash_deg, s_aff)] .= 0
		λ_aff[findall(x->x.p<=trash_deg, λ_aff)] .= 0

		###########################
        # calculate α_aff, μ_aff #
		###########################

        α_aff_pri  = alpha_max(x,x_aff,1.0)
        α_aff_dual = alpha_max(s,s_aff,1.0)
				
		μ = -mean(rxs)
        μ_aff = dot(denoise(x+α_aff_pri*x_aff, tol), denoise(s+α_aff_dual*s_aff, tol))/n
        (μ==0) ? σ = 0 : σ = (μ_aff/μ)^3 # σ = (μ_aff/μ)^3 # 
		
		if show_more
			print("x_aff: "); println(x_aff)
			println("")
			print("s_aff: "); println(s_aff)
			println("")
			print("λ_aff: "); println(λ_aff)
			println("")
			println("")
			print("α_aff_pri: "); println(α_aff_pri)
			print("α_aff_dual: "); println(α_aff_dual)
			println("")
			print("μ_aff: "); println(μ_aff)
			println("")
		end
			
		##############
        # solve 10.7 #
		##############

        rb = zeros(m)
        rc = zeros(n)
        rxs = σ*μ.-α_aff_pri*α_aff_dual*x_aff.*s_aff

        λ_cc,x_cc,s_cc = solve3(f3,rb,rc,rxs)

		x_cc = denoise(x_cc, tol) 
		s_cc = denoise(s_cc, tol)
		λ_cc = denoise(λ_cc, tol)
		x_cc[findall(x->x.p<=trash_deg, x_cc)] .= 0
		s_cc[findall(x->x.p<=trash_deg, s_cc)] .= 0
		λ_cc[findall(x->x.p<=trash_deg, λ_cc)] .= 0
		
		if show_more
			print("x_cc: "); println(x_cc)
			println("")
			print("s_cc: "); println(s_cc)
			println("")
			print("λ_cc: "); println(λ_cc)
			println("")
		end
		
		##############################
        # compute direction and step #
		##############################

        dx = x_aff+x_cc
        dλ = λ_aff+λ_cc
        ds = s_aff+s_cc
		
        α_pri = min(0.99*alpha_max(x,dx,Inf),1)
        α_dual = min(0.99*alpha_max(s,ds,Inf),1)
		
		if show
			print("dx: "); println(dx)
			println("")
			print("ds: "); println(ds)
			println("")
			print("dλ: "); println(dλ)
			println("")
			println("")
			print("α_pri: "); println(α_pri)
			print("α_dual: "); println(α_dual)
			println("")
		end
		
		###############################
        # compute x^k+1, λ^k+1, s^k+1 #
		###############################

        x = x+α_pri*dx
        λ = λ+α_dual*dλ
        s = s+α_dual*ds
		
		x = denoise(x, tol)
		x[findall(x->x<0, x)] .*= -1 #.= 0 #
		
		s = denoise(s, tol)
		s[findall(x->x<0, s)] .*= -1 #.= 0 #
			
		λ = denoise(λ, tol)
		# possible adjusting of λ can be thought (low improvement impact)
		
		# possible cleaning with retrieve_infinitesimals of x, λ and s
		
		###############
        # termination #
		###############

		cost_fun = dot(c,x)+0.5*x'*Q*x
		print("cost: "); println(cost_fun)
		println("")
		r1 = norm(denoise(A*x-b, tol))/rb_den #(1+norm(b))
		r2 = norm(denoise(A'*λ+s-c-Q*x, tol))/rc_den #(1+norm(Q*x+c)) #(1+norm(c)) #
		r3 = denoise(dot(x,s)/n, tol)/(1+abs(cost_fun)) # (magnitude(cost_fun)+abs(cost_fun))
		
		r1 -= retrieve_infinitesimals(r1, trash_deg_r1)
		r2 -= retrieve_infinitesimals(r2, trash_deg_r2)
		r3 -= retrieve_infinitesimals(r3, trash_deg)

        if genLatex
			print("\t$(iter) & \$"); print_latex(mean(x.*s)); print("\$ & \$"); print_latex(x[var_to_show]); print("\$ & \$"); print_latex(cost_fun); println("\$ \\\\");
			println("\t\\hline");
			r = [r
				 r1 r2 r3];
		elseif verbose
            print(iter); print(" "); print(mean(x.*s)); print(" "); print(norm([A'*λ + s - c - Q*x; A*x - b; x.*s])/norm([b;c])); print(" "); print(α_pri); print(" "); println(α_dual); 
        end
		
		if show
			println("")
			print("x: "); println(x)
			println("")
			print("x_std: "); println(map(x->standard_part(x), x))
			println("")
			print("s: "); println(s)
			println("")
			print("λ: "); println(λ)
		end
		
		if show
			println("")
			print("r1: "); println(r1)
			print("r2: "); println(r2)
			print("r3: "); println(r3)
			println("");
			print("μ: "); println(mean(x.*s));
			println("")
		end


        if (typeof(r1)<:Real) ? r1 < tol : all(z->abs(z) < tol, r1.num) 
		
            #r2 = norm(A'*λ+s-c)/(1+norm(c))			

            if (typeof(r2)<:Real) ? r2 < tol : all(z->abs(z) < tol, r2.num) 

                #cx = dot(c,x)
                #r3 = abs(cx-dot(b,λ))/(1+abs(cx))

                if (typeof(r3)<:Real) ? r3 < tol : all(z->abs(z) < tol, r3.num) 

					flag = true;
					if show
						println("")
						print("OPTIMAL SOLUTION X: "); println(x)
						print("OPTIMAL SOLUTION S: "); println(s)
						print("OPTIMAL SOLUTION λ: "); println(λ)
						println("")
					end

					return x,λ,s,flag,iter,r
                end
            end
        end
    end

    return x,λ,s,flag,iter,r
end
