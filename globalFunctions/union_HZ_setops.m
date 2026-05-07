function out = union_HZ_setops(Zvec, options)
% Essentially the "union of many HZs" identity presented in Jake's thesis,
% but with reduced complexity from the union operation with a point
% obtained by simplifying the identity in the CDC/LCSS submission. 
%
% If opt=='sharp', then this returns the sharp union; otherwise, it returns
% a less complex form of the union

n = Zvec(1).n;
for i = 1:length(Zvec)
    if ~isa(Zvec(i), 'abstractZono')
        error(['The input Zvec(' num2str(i) ') must be a zono, conZono, or hybZono object.'])
    end
    
    if i ~= 1 && Zvec(i).n ~= n
        error(['The dimension of Z(' num2str(i) ') does not match dimension of Z(1).'])
    end
    
    if isa(Zvec(i), 'conZono') | isa(Zvec(i), 'zono')
        Z(i) = hybZono(Zvec(i));
    else
        Z(i) = Zvec(i);
    end
end
if exist('options')~=1
    options = '';
end

for i = 1:length(Z)
    Ti(i) = [eye(n); zeros(1, n)]*Z(i) + [zeros(n, 1); 1];
    Ui(i) = union_HZ_vec(Ti(i), zeros(n+1,1), options);

    if i == 1
        SUM = Ui(i);
    else
        SUM = SUM + Ui(i);
    end
end

out = [eye(n) zeros(n,1)]*and(SUM, hybZono(1), [zeros(1,n) 1]);


end
