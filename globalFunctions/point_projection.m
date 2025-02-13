%   Returns the point within zonotopic set Z closest (2-norm) to the point y.
%   Utilizes Gurobi to solve MIQPs/QPs. 
%   
%   Syntax: 
%       [x, xi, result] = point_projection(Z, y)
%   Inputs:
%       Z   - hybZono, conZono, or zono object
%       y   - vector of dimension [Z.n 1]
%   Outputs:
%       x   - point in Z closest to y
%       xi  - factors (xi_c, xi_b) that generate x in Z
%       result - output from GUROBI optimization (for debugging/analysis)
function [x, xi, result] = point_projection(Z, y)    
    if size(y,1)~=Z.n || size(y,2)~=1
        error('Incompatible dimensions; y must be a column vector whose height is the dimensionality of Z.')
    end
    switch class(Z)
        case 'hybZono'
            [x, xi, result] = point_projection_hybZono(Z,y);
        case {'conZono', 'zono'}
            [x, xi, result] = point_projection_zono(Z,y);
    end
end

% implementation for hybZono
function [x, xi, result] = point_projection_hybZono(Z, y)
    n = Z.n; nGc = Z.nGc; nGb = Z.nGb;
    nX = n + nGc + nGb;

    R = eye(n,nX);
    Q = eye(n);
    
    A = zeros(n + Z.nC, nX);
    b = nan(n + Z.nC, 1);
    A(1:n, 1:n) = -eye(n);
    A(1:n, (n+1):(n+nGc)) = Z.Gc;
    A(1:n, (n+nGc+1):(n+nGc+nGb)) = 2*Z.Gb;
    b(1:n) = Z.Gb*ones(nGb,1) - Z.c;

    A((n+1):end, (n+1):(n+nGc)) = Z.Ac;
    A((n+1):end, (n+nGc+1):(n+nGc+nGb)) = 2*Z.Ab;
    b((n+1):end) = Z.b + Z.Ab*ones(nGb,1);

    z_lb = zeros(n,1);
    z_ub = zeros(n,1);
    basisVectors = eye(n);
    for i = 1:n
        [s,~] = supportFunc(Z,-basisVectors(:,i));
        z_lb(i) = -s;
        [s,~] = supportFunc(Z,basisVectors(:,i));
        z_ub(i) = s;
    end

    lb = [z_lb; -ones(nGc,1); zeros(nGb, 1)];
    ub = [z_ub; ones(nGc+nGb,1)];

    vType(1:(n+nGc)) = 'C';
    vType((n+nGc+1):(n+nGc+nGb)) = 'B';
    
    model.A = sparse(A);
    model.obj = -2*y'*Q*R;
    model.sense = '=';
    model.rhs = b;
    model.lb = lb;
    model.ub = ub;    
    model.vtype = vType;
    model.modelsense = 'min';
    model.objcon = y'*Q*y;
    model.Q = sparse(R'*Q*R);

    % GUROBI settings
    params.Threads = 8;
    params.outputFlag = 0;
    params.PoolSearchMode = 1;
    params.PoolSolutions = 2e9;
    params.MIPFocus = 1;
    result = gurobi(model, params);

    x = result.x(1:n);
    xi = result.x(n+1:end);
    xi((nGc+1):end) = 2*xi(nGc+1:end)-1;
end

% used in prototyping
function [x,xi,result] = point_projection_hybZono_yalmip(Z,y)
    n = Z.n; nGc = Z.nGc; nGb = Z.nGb;

    x = sdpvar(n,1);
    xi_c = sdpvar(nGc,1);
    xi_b = binvar(nGb,1);
    opts = sdpsettings;
    opts.solver = 'gurobi';
    opts.gurobi.MIPFocus = 0;
    opts.gurobi.OutputFlag = 0;
    opts.gurobi.Threads = 8;
    opts.gurobi.OptimalityTol = 1e-8;
    opts.verbose = 0;
    J = norm(x-y,2) + 1e-4*norm([xi_c; xi_b], 2);
    cons = [ norm(xi_c, 'inf') <= 1 ,
             x == Z.Gc*xi_c + Z.Gb*(2*xi_b - ones(nGb,1)) + Z.c,
             Z.b == Z.Ac*xi_c + Z.Ab*(2*xi_b - ones(nGb,1))];
    result = optimize(cons,J,opts);
    x = value(x);
    xi = [value(xi_c); 2*value(xi_b)-1];
end

function [x, xi, result] = point_projection_zono(Z,y)
    n = Z.n; nG = Z.nG;
    nX = n + nG;

    R = eye(n,nX);
    Q = eye(n);
    
        A = zeros(n + Z.nC, nX);
        b = nan(n + Z.nC, 1);
        A(1:n, 1:n) = -eye(n);
        A(1:n, (n+1):end) = Z.G;
        b(1:n) = -Z.c;
    if isa(Z, 'conZono')
        A((n+1):end, (n+1):(end)) = Z.A;
        b((n+1):end) = Z.b;
    end
        
    z_lb = zeros(n,1);
    z_ub = zeros(n,1);
    basisVectors = eye(n);
    for i = 1:n
        [s,~] = supportFunc(Z,-basisVectors(:,i));
        z_lb(i) = -s;
        [s,~] = supportFunc(Z,basisVectors(:,i));
        z_ub(i) = s;
    end

    lb = [z_lb; -ones(nG,1)];
    ub = [z_ub; ones(nG,1)];

    vType(1:(n+nG)) = 'C';
    model.A = sparse(A);
    model.obj = -2*y'*Q*R;
    model.sense = '=';
    model.rhs = b;
    model.lb = lb;
    model.ub = ub;    
    model.vtype = vType;
    model.modelsense = 'min';
    model.objcon = y'*Q*y;
    model.Q = sparse(R'*Q*R);

    params.Threads = 8;
    params.outputFlag = 0;
    params.PoolSearchMode = 1;
    params.PoolSolutions = 2e9;
    params.MIPFocus = 0;
    result = gurobi(model, params);

    x = result.x(1:n);
    xi = result.x(n+1:end);
end