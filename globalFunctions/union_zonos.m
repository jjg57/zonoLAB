function Z = union_zonos(Z_array)
% union_zonos takes an array of zonotopes (strictly zonotopes, no
% constraints) and efficiently represents their union as a hybrid zonotope.
% All of the generator matrices are checked to find the unique generators
% and remove duplicates between sets. This union is sharp.
% 
% Syntax: 
%   Z = union_zonos(Z_array)
%
% Inputs:
%   Z_array - a vector array of zonos [Z1 Z2 ... ZN]
%
% Outputs:
%   Z - a hybZono that represents the sharp union of Z1, Z2, ..., ZN
% 

% to-do: input checking
N = length(Z_array);
n = Z_array(1).n;
nG_arr = [];
G = [];
C = [];
for i = 1:N
    if Z_array(i).n~=n
        error('union_zono:inconsistent_dims', 'Dimensions of sets in array are inconsistent.')
    end
    Z_array(i) = to_0_1_form(Z_array(i));
    nG_arr(i) = Z_array(i).nG; % store number of generators in Z_i
    G = [G Z_array(i).G]; % store all generators
    C = [C Z_array(i).c]; % store all centers
end

M = []; % incidence matrix

[Gtilde, unique_id, col_id] = unique(G', 'rows', 'stable'); % find the unique generators
Gtilde = Gtilde';
nG = length(unique_id); % how many unique gens are there

for i = 1:nG
    Ntilde_i(i) = sum(col_id==i); % count how many times each unique generator appears
end
Ntilde = diag(Ntilde_i);

idx = 0;
for i = 1:N
    this_col_id = col_id(idx+1 : idx+nG_arr(i));
    m = zeros(nG,1);
    m(this_col_id) = 1; % make this column of the incidence matrix
    M = [M m];

    idx = idx + nG_arr(i);
end

Z = hybZono([Gtilde zeros(n,nG)], C, zeros(n,1), [eye(nG) Ntilde; zeros(1,2*nG)], [-M; ones(1,N)], [zeros(nG,1); 1]);

Z = to_m1_1_form(Z);
end
