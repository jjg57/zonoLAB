%	interval2hybZono - convert an interval set to hybzono
% 
%	Syntax: 
%		Zh = interval2hybZono(box)
% 
%	Inputs:
%		box : - n x 2 dimensional matrix containing bounds on each of the 
%				dimensions -> box(:,1) = lb , box(:,2) = ub
%			  - Polyhedron MPT3 object, if this poly is not an interval 
%				then this function will return its intervalbox in HCG-rep
% 
%	Outputs:
%		Zh : n dimensional hybrid zonotope in HCG-rep
%	
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
%   Modified for use with zonoLAB by Jonah Glunt - jglunt@psu.edu

function Zh = interval2hybZono(box)

% if box is a polyhedron, convert it to an interval
if isa(box,'Polyhedron')
	Pbox = box;
	Pbox = Pbox.outerApprox();
	box = [];
	box(:,1) = Pbox.Internal.lb;
	box(:,2) = Pbox.Internal.ub;
elseif size(box,2) > 2
	error('Input must be a MPT3 polyhedron or an n x 2 dimensional matrix containing bounds on each of the dimensions -> box(:,1) = lb , box(:,2) = ub')
end

% intervals are simple zonotopes, can represent them in HCG-rep
n = size(box,1);
Gc = diag((box(:,2)-box(:,1))/2);
c = (box(:,2)+box(:,1))/2;
% Zh = hybZono(Gc,zeros(n,0),c,zeros(n,0),zeros(n,0),zeros(n,0));
Zh = hybZono(Gc,[],c,[],[],[]);
% and remove zeros if some dimensions are a single point
try
    Zh = removeZeros(Zh);
catch
    warning('The function "removeZeros(Zh)" is not yet ported for zonoLAB compatability.')
end

end