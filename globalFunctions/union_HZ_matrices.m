function Zh = union_HZ_matrices(HZ_vec, options)
% The union method presented in initial submission to CDC/L-CSS. The
% operation is defined by constructing the matrices (as opposed to the
% other version, which uses a simpler sharp union of HZ and a vector and
% other sharpness preserving set operations such as Mink sum and hyperplane
% intersection). 
%
% If opt=='sharp', then this returns the sharp union; otherwise, it returns
% a less complex form of the union
    
    % convert all sets to hybzonos if they are not already
    for ii = 1:length(HZ_vec)
        if isa(HZ_vec(ii), 'zono')
            HZ_vec_new(ii) = hybZono(HZ_vec(ii).G, zeros(HZ_vec(ii).n, 0), HZ_vec(ii).c, zeros(0, HZ_vec(ii).nG), zeros(0, 0), zeros(0,1));
        elseif isa(HZ_vec(ii), 'conZono')
            HZ_vec_new(ii) = hybZono(HZ_vec(ii).G, zeros(HZ_vec(ii).n, 0), HZ_vec(ii).c, HZ_vec(ii).A, zeros(HZ_vec(ii).nC, 0), HZ_vec(ii).b);
        elseif ~isa(HZ_vec(ii), 'hybZono')
            error('unknown data type')
        else
            HZ_vec_new(ii) = HZ_vec(ii);
        end
    end
    HZ_vec = HZ_vec_new;

    % need generators in [0,1] range
    for ii = 1:length(HZ_vec)
        HZ_vec(ii) = to_0_1_form(HZ_vec(ii));
    end
    
    if exist('options')~=1
        options = '';
    end
    %% build hybzono
    n = HZ_vec(1).n;
    num_zonos = length(HZ_vec);

    if isequal(options, 'sharp')

        % generators
        Gc = zeros(n,0); % init
        Gb = zeros(n, 0); % init
        nG = nan(1,num_zonos);
        for ii = 1:num_zonos
            nG(ii) = HZ_vec(ii).nGc+HZ_vec(ii).nGb;
            Gc = [Gc, HZ_vec(ii).Gc, zeros(n,nG(ii))];
            Gb = [Gb, HZ_vec(ii).Gb, HZ_vec(ii).c];
        end
    
        c = zeros(n, 1);
        
        % constraints
        Ac = []; % init
        Ab = []; % init
        % b = []; % init
        ind_b = [];
        idx = 1;
        Ab_f = [];
        for ii = 1:num_zonos
            Ac_ii = [HZ_vec(ii).Ac , zeros(HZ_vec(ii).nC, nG(ii));
                     [eye(HZ_vec(ii).nGc) ; zeros(HZ_vec(ii).nGb, HZ_vec(ii).nGc)] , eye(nG(ii))];
            Ab_ii = [HZ_vec(ii).Ab , -HZ_vec(ii).b ; 
                     [zeros(HZ_vec(ii).nGc, HZ_vec(ii).nGb) ; eye(HZ_vec(ii).nGb)] , -ones(nG(ii),1)];
            % b_ii = [0; zeros(HZ_vec(ii).nC,1)];
       
            if isempty(Ac)
                Ac = Ac_ii;
                Ab = Ab_ii;
            else
                Ac = blkdiag(Ac, Ac_ii);
                Ab = blkdiag(Ab, Ab_ii);
            end
    
            Ab_f = [Ab_f zeros(1,HZ_vec(ii).nGb) 1];
        end
        
        % choose one
        Ac = [Ac ; zeros(1,size(Ac,2))];
        Ab = [Ab ; Ab_f];
        b = zeros(size(Ac,1),1);
        b(end) = 1;
    else
        Gc = zeros(n,0); % init
        for ii = 1:num_zonos
            Gc = [Gc, HZ_vec(ii).Gc, zeros(n,1)];
        end

        Gb = zeros(n, 0); % init
        for ii = 1:num_zonos
            Gb = [Gb, HZ_vec(ii).Gb, HZ_vec(ii).c];
        end

        c = zeros(n, 1);

        % constraints
        Ac = []; % init
        Ab = []; % init
        b = []; % init
        ind_b = [];
        idx = 1;
        for ii = 1:num_zonos
            Ac_ii = [ones(1,HZ_vec(ii).nGc), HZ_vec(ii).nGc+HZ_vec(ii).nGb;
                     HZ_vec(ii).Ac, zeros(HZ_vec(ii).nC,1)];
            Ab_ii = [ones(1,HZ_vec(ii).nGb), -HZ_vec(ii).nGc-HZ_vec(ii).nGb;
                     HZ_vec(ii).Ab, -HZ_vec(ii).b];
            b_ii = [0; zeros(HZ_vec(ii).nC,1)];

            if isempty(Ac)
                Ac = Ac_ii;
                Ab = Ab_ii;
                b = b_ii;
            else
                Ac = [Ac, zeros(size(Ac,1), size(Ac_ii,2));
                      zeros(size(Ac_ii,1), size(Ac,2)), Ac_ii];
                Ab = [Ab, zeros(size(Ab,1), size(Ab_ii,2));
                      zeros(size(Ab_ii,1), size(Ab,2)), Ab_ii];
                b = [b;
                    b_ii];
            end

            % binary indices tracking
            ind_b = [ind_b, idx + HZ_vec(ii).nGb];
            idx = idx + HZ_vec(ii).nGb + 1;
        end

        % choose one
        Ac = [Ac;
            zeros(1,size(Ac,2))];
        Ab = [Ab;
            zeros(1,size(Ab,2))];
        Ab(end,ind_b) = 1;

        b = [b;
            1];
    end
    
    % convert to [-1,1] form
    Zh = hybZono(Gc, Gb, c, Ac, Ab, b);
    Zh = to_m1_1_form(Zh);
   
    
    end