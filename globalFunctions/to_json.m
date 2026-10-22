function to_json(Z, filename)
%TO_JSON Write a hybrid zonotope to a JSON file.
    
    Z = simplifySetType(Z);
    data = struct();

    %% Determine class
    if isa(Z, 'hybZono')
        data.class = "HybZono";
    elseif isa(Z, 'conZono')
        data.class = "ConZono";
        Z = hybZono(Z.G, [], Z.c, Z.A, [], Z.b);
    elseif isa(Z, 'zono')
        data.class = "Zono";
        Z = hybZono(Z.G, [], Z.c, [], [], []);
    else
        error("to_json: Unrecognized set type.");
    end

    %% Serialize

    data.n  = Z.n;

    data.Gc = sparse_to_json(Z.Gc);
    data.Gb = sparse_to_json(Z.Gb);

    data.c = Z.c();

    data.Ac = sparse_to_json(Z.Ac);
    data.Ab = sparse_to_json(Z.Ab);

    data.b = Z.b;

    data.zero_one_form = false;

    %% Write file

    fid = fopen(filename,'w');

    if fid==-1
        error("Could not open file for writing.");
    end

    fwrite(fid,jsonencode(data),'char');
    fclose(fid);

end