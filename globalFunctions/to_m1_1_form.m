function Z_m11 = to_m1_1_form(Z)

if isa(Z, 'zono')
    c = Z.c + 0.5*Z.G*ones(Z.nG, 1);
    G = 0.5*Z.G;
    Z_m11 = zono(G,c);
elseif isa(Z, 'conZono')
    c = Z.c + 0.5*Z.G*ones(Z.nG, 1);
    G = 0.5*Z.G;
    A = 0.5*Z.A;
    b = Z.b - 0.5*Z.A*ones(Z.nG, 1);
    Z_m11 = conZono(G,c,A,b);
elseif isa(Z, 'hybZono')
    c = Z.c + 0.5*[Z.Gc, Z.Gb]*ones(Z.nGc + Z.nGb, 1);
    b = Z.b - 0.5*[Z.Ac, Z.Ab]*ones(Z.nGc + Z.nGb, 1);
    Gc = 0.5*Z.Gc;
    Gb = 0.5*Z.Gb;
    Ac = 0.5*Z.Ac;
    Ab = 0.5*Z.Ab;
    Z_m11 = hybZono(Gc, Gb, c, Ac, Ab, b);
else
    error('invalid data type for Z')
end

end