function Zout = mpower(Z, n)
    % Cartesian product of Z with itself n times: Z^n
    if ~isnumeric(n) || ~isscalar(n) || n < 1 || n ~= floor(n)
        error('Exponent must be a positive integer.');
    end

    Zout = Z;
    for i = 2:n
        Zout = cartProd(Zout, Z);
    end
end