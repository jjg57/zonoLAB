%	interval2zono - convert interval bounds to a zonotope
%   Can be used instead of the existing function, interval2hybZono, which
%   returns a hybZono object whose entries for Gb, Ac, Ab, and b. This
%   function instead returns a pure zono object. 
% 
%	Syntax: 
%		Zh = interval2hybZono(box)
% 
%	Inputs:
%		box : - n x 2 dimensional matrix containing bounds on each of the 
%				dimensions -> box(:,1) = lb , box(:,2) = ub
% 
%	Outputs:
%		Z : n dimensional hybrid zonotope in G-rep
%	
%	Jonah Glunt - jglunt@psu.edu

function Zh = interval2zono(box)

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
Zh = zono(Gc,c);
end