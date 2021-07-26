include("../src/ipqp.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/createTable.jl")

# TODO implement test suit to verify correctness of NA-IPM implementation on all the benchmarks

function add_slack(c,Q, A)
	A = [A, I]
	sA1 = size(A,1)
	sQ1 = size(Q,1)
	c = vcat(c, zeros(sA1))
	Q = [Q zeros(sQ1, sA1); zeros(sA1, sQ1+sA1)]

	return c, Q, A
end

# linear lex kite

function llk(z)
	Q = zeros(2,2)

	if z
		c = [-8-14η, -12-10η]; 	# converges to (30, 50)
	else
		c = [-8-4η, -12-10η]; 	# converges to (0, 70)
	end
	
	b = [120, 210, 270, 60];

	A = [ 2  1;
		  2  3;
		  4  3;
		 -1 -2];

	c, Q, A = add_slack(c, Q, A)

	return c, Q, A, b
end

# quadratic lex kite

function qlk(z)

	if z
		Q = [1 0; 0 1]*η # converges to (30, 50)
	else
		Q = [-1 0; 0 -1]*η # converges to (0, 70)
	end
	
	c = [-8, -12]
	
	b = [120, 210, 270, 60];

	A = [ 2  1;
		  2  3;
		  4  3;
		 -1 -2];
		 
	c, Q, A = add_slack(c, Q, A)
	
	return c, Q, A, b
end

# linear NA kite

function lnk()

	Q = zeros(2,2)
	
	c = [-3-η, -1+2η]
	
	A = [ 2     1-η;
		  2+η   3+2η;
		  4+2η  3-η;
		 -1+η  -2+3η];
	
	b = [120-3η, 210-η, 270-4η, -60-η];
	
	c, Q, A = add_slack(c, Q, A)
	
	return c, Q, A, b
end

# quadratic NA kite

function qnk()

	Q = [2+η -3-2η; -3-2η 10-3η].*η;
	
	c = [-3-η, -1+2η]
	
	A = [ 2     1-η;
		  2+η   3+2η;
		  4+2η  3-η;
		 -1+η  -2+3η];
	
	b = [120-3η, 210-η, 270-4η, -60-η];
	
	c, Q, A = add_slack(c, Q, A)
	
	return c, Q, A, b
end

# 2obj lex pyramid

function 2lp()

	Q = [10 -2  4;
		 -2 10  4;
		  4  4  4];

	c = [-16-η, -16-η, -16, 0, 0, 0, 0]; # converges to (1.5,1.5,0)

	b = [1, -1, 1, 3];

	A = [-1  1  1;  # the four faces of the pyramid
		 -1 -1  1;
		  1 -1  1;
		  1  1  1];

	c, Q, A = add_slack(c, Q, A)
	
	return c, Q, A, b
end 

# 3obj lex pyramid

function 3lp()

	Q = [10 -2  4;
		 -2 10  4;
		  4  4  4];

	c = [-1-10η-5η*η, -1-10η-3η*η, -1-2η*η, 0, 0, 0, 0]; # converges to (1.67, 1.17, 0.17)

	b = [1, -1, 1, 3];

	A = [-1  1  1;  # the four faces of the pyramid
		 -1 -1  1;
		  1 -1  1;
		  1  1  1];

	c, Q, A = add_slack(c, Q, A)
	
	return c, Q, A, b
end 

# lex distance polyhedra

function embed_2D(A_x, A_y, type="point";x=nothing,y=nothing)
	
	s_A_x = size(A_x)
	s_A_y = size(A_y)
	num_slack = s_A_x[1] + s_A_y[1];
	A = [A_x zeros(s_A_x[1], s_A_y[2]); zeros(s_A_y[1], s_A_x[2]) A_y];
	A = [A I];
	
	if type=="point"
		(x==nothing || y==nothing) && throw(ArgumentError("x and y must be specified in type point"))
		
		c = vcat(-2η*vcat(x,y), zeros(num_slack))
		Q = 2*[1+η  0  -1   0;
		        0  1+η  0  -1;
			   -1   0  1+η  0;
			    0  -1   0  1+η];
		Q = [Q zeros(4,num_slack);
		     zeros(num_slack, 4+num_slack)];
	elseif type=="barycenter"
		throw(ArgumentError("Not implemented yet"))
	end

	return A, Q, c
end

function ldp()
	A_x = [-1  0;
		   -1  1;
			1  1;
			1 -1;
			0 -1];

	A_y = [-3 -1;
			1  3;
			1  1;
			3  1;
		   -1 -3;
		   -1 -1];

	b = [-1, 3, 10, 2, -1.5, -25, 35, 17, 40, -20, -13];

	x_ref = [2,5];
	y_ref = [10.5, 4.5];

	A, Q, c = embed_2D(A_x, A_y, "point", x=x_ref, y=y_ref)
	
	return c, Q, A, b
end

# 2obj lex floor planning

function constraints3(W, H)
	# horizontal positioning
	A_h = [1  0 -1  0  0  0;
		   0  1 -1  0  0  0;
		   0  0  1  0  0  0];
	 
	b_h = [0; 0; W];


	# vertical positioning
	A_v = [0  0  0 -1  1  0;
		   0  0  0  1  0  0;
		   0  0  0  0  0  1];
	 
	b_v = [0; H; H];


	# upper aspect ratio
	A_u = [2  0 -2 -1  0  0;
		   0  2 -2  1 -1  0;
		   0  0  2  0  0 -1];
	   
	b_u = [-H; 0; 2*W-H];


	# lower aspect ratio
	A_l = [-1  0  1  2  0  0;
		    0 -1  1 -2  2  0;
		    0  0 -1  0  0  2];

	b_l = [2*H; 0; 2*H-W];


	# minimal height
	A_m = [0  0  0  1  0  0;
		   0  0  0 -1  1  0;
		   0  0  0  0  0  1];
	   
	b_m = [H-1; -1; H-1];

	A = [A_h; A_v; A_u; A_l; A_m];
	
	b = [b_h; b_v; b_u; b_l; b_m];

	return A, b
end

function cost_f3(W, H)
	Q = -[0  0  0  1  0  0;
	      0  0  0 -1  1  0;
		  0  0  0  0 -1  1;
		  1 -1  0  0  0  0;
		  0  1 -1  0  0  0;
		  0  0  1  0  0  0]
	Q+= [ 4 -2 -4  0  0  0;
	     -2  2  0  0  0  0;
		 -4  0  8  0  0  0;
		  0  0  0 10 -4 -2;
		  0  0  0 -4  2  0;
		  0  0  0 -2  0  2].*η;
	
	c = -[-H, 0, 0, 0, 0, -W];
	c+=  [2*W, 0, -4*W, -4*H, 2*H, 0]*η
	
	return Q, c
end

function 2lfp(H, W)
	H = 10
	W = 14
	
	A, b = constraints3(W,H)

	Q, c = cost_f3(W,H)
	
	c, Q, A = add_slack(c, Q, A)
	
	return A, b, c, Q
end

# ---------------------------------------------- #

function print_result(flag, text)
		print(text);
		ifelse(flag, println(" success"), println(" failed")) 
		println("")
end

tol=1e-8;
genLatex = false;
verbose = false;

println("Test launched, please wait...")

c, Q, A, b = llk(true);
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);
print_result(sol.flag, "llk_t:")

c, Q, A, b = llk(false);
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);
print_result(sol.flag, "llk_f:")

c, Q, A, b = qlk(true);
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);
print_result(sol.flag, "qlk_t:")

c, Q, A, b = qlk(false);
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);
print_result(sol.flag, "qlk_f:")

c, Q, A, b = lnk();
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);
print_result(sol.flag, "lnk:")

c, Q, A, b = qnk();
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=3:6);
print_result(sol.flag, "qnk:")

c, Q, A, b = 2lp();
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=4:8);
print_result(sol.flag, "2np:")

c, Q, A, b = 2lp();
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=4:8);
print_result(sol.flag, "2lp:")

c, Q, A, b = ldp();
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=5:size(A,2));
print_result(sol.flag, "ldp:")

c, Q, A, b = 2lfp(10, 14);
sol = ipqp(A,b,c,Q, tol; maxit=30, verbose=verbose, genLatex=genLatex, slack_var=5:size(A,2));
print_result(sol.flag, "2lfp:")