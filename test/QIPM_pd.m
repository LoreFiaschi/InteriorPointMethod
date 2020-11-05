H = -2*eye(2);
f = [2; 2];
A = [-1 1; 1 1; 1 1; -1 1];
b = [-1; 1; 3; -1];

n = length(f);
Aeq = zeros(0,n);
beq = zeros(0,1);

x0 = [1.22; 1.22];

opt = mpcInteriorPointOptions;

[x,exitflag] = mpcInteriorPointSolver(H,f,A,b,Aeq,beq,x0,opt);

x