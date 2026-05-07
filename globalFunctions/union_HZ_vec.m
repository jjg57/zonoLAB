function out = union_HZ_vec(Z, v, options)
% The efficient (and sharp) union of a hybrid zonotope Z with a given
% vector v
% 
% If opt=='sharp', then this returns the sharp union; otherwise, it returns
% a less complex form of the union

if ~isa(Z, 'abstractZono')
    error('The input Z must be a zono, conZono, or hybZono object.')
end

if ~isa(v, 'double') & isequal(size(v), [Z.n 1])
    error('The dimension of the set and vector do not match.')
end

if isa(Z, 'conZono') | isa(Z, 'zono')
    Z = hybZono(Z);
end
X = to_0_1_form(Z);

if exist('options')~=1
    options = '';
end

if isequal(options, 'sharp')
    Gc = [X.Gc zeros(X.n, X.nGc+X.nGb)];
    Gb = [X.Gb X.c-v];
    c = v;
    Ac = [X.Ac zeros(X.nC,X.nGc+X.nGb) ; 
         [eye(X.nGc); zeros(X.nGb,X.nGc)] eye(X.nGc+X.nGb)];
    Ab = [X.Ab -X.b; 
         [zeros(X.nGc, X.nGb); eye(X.nGb)] -ones(X.nGc+X.nGb,1)];
    b = [zeros(X.nC,1); zeros(X.nGc+X.nGb,1)];
    out = hybZono(Gc, Gb, c, Ac, Ab, b);
else
    Gc = [X.Gc zeros(X.n,1)];
    Gb = [X.Gb, X.c-v];
    c = v;
    Ac = [X.Ac , zeros(X.nC, 1);
          ones(1,X.nGc) , X.nGc+X.nGb];
    Ab = [X.Ab , -X.b;
          ones(1,X.nGb) , -X.nGc-X.nGb];
    b = zeros(X.nC+1,1);
    out = hybZono(Gc, Gb, c, Ac, Ab, b);
end
out = to_m1_1_form(out);
end