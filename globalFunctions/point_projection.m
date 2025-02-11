function x = point_projection(Z, y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    optSolver = solverOptions;
    opSolver.milpSolver = 'gurobi';
    n = Z.n; nGc = Z.nGc; nGb = Z.nGb;
    nX = n + nGc + nGb;

    R = eye(n,nX);
    Q = eye(n);

    Aeq = [Z.Ac 2*Z.Ab];
    beq = Z.b + Z.Ab*ones(Z.nGb,1);
    ub = ones(Z.nGc+Z.nGb, 1);
    lb = [ones(Z.nGc, 1); zeros(Z.nGb, 1)];
    vType(1:Z.nGc) = 'C';
    vType(Z.nGc+1:(Z.nGc+Z.nGb)) = 'B';

    model.sense = '=';
    
end