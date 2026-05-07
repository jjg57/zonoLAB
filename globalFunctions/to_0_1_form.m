function Z_01 = to_0_1_form(Z)

if isa(Z, 'zono')
    c = Z.c - Z.G*ones(Z.nG, 1);
    G = 2*Z.G;
    Z_01 = zono(G,c);
elseif isa(Z, 'conZono')
    c = Z.c - Z.G*ones(Z.nG, 1);
    G = 2*Z.G;
    A = 2*Z.A;
    b = Z.b + Z.A*ones(Z.nG, 1);
    Z_01 = conZono(G,c,A,b);
elseif isa(Z, 'hybZono')
    c = Z.c - [Z.Gc, Z.Gb]*ones(Z.nGc + Z.nGb, 1);
    b = Z.b + [Z.Ac, Z.Ab]*ones(Z.nGc + Z.nGb, 1);
    Gc = 2*Z.Gc;
    Gb = 2*Z.Gb;
    Ac = 2*Z.Ac;
    Ab = 2*Z.Ab;
    Z_01 = hybZono(Gc, Gb, c, Ac, Ab, b);
else
    error('invalid data type for Z')
end

end