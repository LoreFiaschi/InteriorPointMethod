using LinearAlgebra
M = 30;

A = [zeros(1,4) one(Float64);
	 Matrix{Float64}(I,4,4) M.*[ones(2,1); -ones(2,1)];
	-Matrix{Float64}(I,2,2) zeros( 2, 2) M.*ones(2,1)];

A_dom = [zeros( 1, 2) -1  0  2;
		 zeros( 1, 2)  0 -1  4;
		 -1  1 zeros( 1, 2)  1];

A = [A; A_dom];

A