function result = point_projection(Z, y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    optSolver = solverOptions;
    opSolver.milpSolver = 'gurobi';
    n = Z.n; nGc = Z.nGc; nGb = Z.nGb;
    nX = n + nGc + nGb;

    R = eye(n,nX);
    Q = eye(n);
    
    A = zeros(n + Z.nC, nX);
    b = nan(n + Z.nC, 1);
    A(1:n, 1:n) = -eye(n);
    A(1:n, (n+1):(n+nGc)) = Z.Gc;
    A(1:n, (n+nGc+1):(n+nGc+nGb)) = 2*Z.Gb;
    b(1:n) = -Z.Ab*ones(nGb,1);

    A((n+1):end, (n+1):(n+nGc)) = Z.Ac;
    A((n+1):end, (n+nGc+1):(n+nGc+nGb)) = 2*Z.Ab;
    b((n+1):end) = Z.b - Z.Ab*ones(nGb,1);

    z_lb = zeros(obj.n,1);
    z_ub = zeros(obj.n,1);
    basisVectors = eye(obj.n);
    for i = 1:obj.n
        [s,~] = supportFunc(obj,-basisVectors(:,i));
        z_lb(i) = -s;
        [s,~] = supportFunc(obj,basisVectors(:,i));
        z_ub(i) = s;
    end

    lb = [z_lb; -ones(nGc,1); zeros(nGb, 1)];
    ub = [z_ub; ones(nGc+nGb,1)];

    vType(1:Z.nGc) = 'C';
    vType(Z.nGc+1:(Z.nGc+Z.nGb)) = 'B';
    
    model.A = A;
    model.obj = -2*y'*Q*R;
    model.sense = '=';
    model.rhs = b;
    model.lb = lb;
    model.ub = ub;    
    model.vtype = vType;
    model.modelsense = 'min';
    model.objcon = y'*Q*y;
    model.Q = R'*Q*R;

    result = gurobi(model, params);
end