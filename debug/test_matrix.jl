using LinearAlgebra
M = 30;

A = [zeros(Float64,1,4) one(Float64);
	 Matrix{Float64}(I,4,4) M.*[ones(Float64,2,1); -ones(Float64,2,1)];
	-Matrix{Float64}(I,2,2) zeros(Float64, 2, 2) M.*ones(Float64,2,1)];

A_dom = [-1  1 zeros(Float64, 1, 2) -1;
		  1 -1 zeros(Float64, 1, 2) -1;
		  zeros(Float64, 2, 2)  Matrix{Float64}(I,2,2)    [-2; -2];
		  zeros(Float64, 2, 2) -Matrix{Float64}(I,2,2)    [ 1;  1] ];

A = [A; A_dom];

A
