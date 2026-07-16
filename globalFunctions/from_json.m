function Z = from_json(filename)
%FROM_JSON Read a hybrid zonotope from a JSON file.

    fid = fopen(filename,'r');

    if fid==-1
        error("Could not open file.");
    end

    raw = fread(fid,'*char')';
    fclose(fid);

    data = jsondecode(raw);

    %% Deserialize

    Gc = json_to_sparse(data.Gc);
    Gb = json_to_sparse(data.Gb);

    c = data.c(:);

    Ac = json_to_sparse(data.Ac);
    Ab = json_to_sparse(data.Ab);

    b = data.b(:);

    zero_one_form = data.zero_one_form;

    n = data.n;

    %% Construct object

    switch char(data.class)

        case 'HybZono'
            Z = hybZono(Gc,Gb,c,Ac,Ab,b);

        case 'ConZono'
            Z = conZono(Gc,c,Ac,b);

        case 'Zono'
            Z = zono(Gc,c);

        otherwise
            error("Unknown class '%s'.",data.class);

    end

    if zero_one_form
        Z = to_m1_1_form(Z);
    end

end