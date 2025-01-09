%% First runs the code from the Logical Zonotope Toolbox, Then the HZ comparison

clear all

rng(100);
steps = 10;
dim =20;
numOfVar = 20;
numOfVarU = 20;
numOfPoints = 10;
for i=1:numOfVar
    Apt{i} = randi([0 1],dim,numOfPoints);
    A{i} = reduce(logicalZonotope.enclosePoints(Apt{i}));
end

for i=1:numOfVarU
    pointsU{i} = randi([0 1],dim,numOfPoints);
    U{i} = reduce(logicalZonotope.enclosePoints(pointsU{i}));
end
x=1;
tic
for i =1:steps
    [A] = boolFun(A,U);
end
execlogicZono = toc

%% ----------- Hybrid Zonotopes ---------

plot_LOGIC = 1;
Add4Plot = hybZono(zeros(3,1),.1*eye(3),[]);


%---------------- XNOR (not most efficient but what-evs) ---------------
XNOR = hybZono([.5;.5;0],zeros(3,1),[.5;-.5;0]);
XNOR = union(XNOR,hybZono([.5;.5;1],[],[.5;.5;0]));

XNORp = XNOR + Add4Plot;

if plot_LOGIC
figure
subplot(221)
grid on; grid MINOR
XNORp.plot('method','MPT')
title('XNOR')
xlabel('x')
ylabel('y')
end

% ----------------------- AND (not most efficient but what-evs) ----------
AND = hybZono([0;.5;0],zeros(3,1),[0;-.5;0]);
AND = union(AND,hybZono([1;.5;.5],[],[0;.5;.5]));

ANDp = AND + Add4Plot;

if plot_LOGIC
% figure
subplot(222)
grid on; grid MINOR
ANDp.plot('method','MPT')
title('AND')
xlabel('x')
ylabel('y')
end

% ----------- NAND ------------------------------------------------------- 
NAND = diag([1 1 -1])*AND + [0 0 1]';

NANDp = NAND + Add4Plot;

if plot_LOGIC
% figure
subplot(223)
grid on; grid MINOR
NANDp.plot('method','MPT')
title('NAND')
xlabel('x')
ylabel('y')
end

%--------------------- OR (not most efficient but what-evs)----------------
OR = hybZono([1;.5;1],zeros(3,1),[0;.5;0]);
OR = union(OR,hybZono([0;.5;.5],[],[0;.5;.5]));

ORp = OR + Add4Plot;

if plot_LOGIC
% figure
subplot(224)
grid on; grid MINOR
ORp.plot('method','MPT')
title('OR')
xlabel('x')
ylabel('y')
end

%---------------------- Let us make some boolean functions :) ------------
% See the boolean functions from the "High-Dimensional Boolean Function"
% in the paper

% Create basis
Phi = hybZono(.5*ones(12,1),zeros(12,1),.5*eye(12));
% Couple the intermediate variables
% [b1, b2, b3, u1, u2, u3, tmp1=xnor(b2,b1), tmp2=and(b1,u2), tmp3=xnor(u2,u3), b1+, b2+, b3+]^T
Phi = andR(Phi,XNOR,M_i([2 1 7],12));
Phi = andR(Phi,AND ,M_i([1 5 8],12));
Phi = andR(Phi,XNOR,M_i([5 6 9],12));
% Couple k+1 variables
% [b1, b2, b3, u1, u2, u3, tmp1, tmp2, tmp3, b1+=or(u1,tmp1), b2+=xnor(b2,tmp2), b3+=nand(b3,tmp3)]^T
Phi = andR(Phi,  OR,M_i([4 7 10],12));
Phi = andR(Phi,XNOR,M_i([2 8 11],12));
Phi = andR(Phi,NAND,M_i([3 9 12],12));
% Tranform to get rid of intermediate states
Phi = M_i([1:6, 10:12],12)*Phi

%-------------- Test Phi and Compare Results to the Boolean Function (See subfunctions)

m = 0;
for i = 1:100
R3 = round(rand(3,1));
R3set = hybZono(R3,[],[]);
U = round(rand(3,1));
Uset = hybZono(U,[],[]);

% Successor Set Identity
R3plusSet = [zeros(3,6) eye(3)]*andR(Phi,cartProd(R3set,Uset),[eye(6) zeros(6,3)]);

% Note Gc=0, and # leaves = 1
R3plus = R3plusSet.c+R3plusSet.Gb*R3plusSet.getLeaves;

R3 = logical(R3);
U = logical(U);

% UNCOMMENT THE FOLLOWING LINE IF YOU WANT TO VERIFY ACCURACY BY HAND
[R3, U, R3plus, BooleanFunction(R3,U)];

m = m + max(abs(R3plus-BooleanFunction(R3,U)));
disp('Testing State-Update Set Logic 100x')
disp('[Iteration, Total Deviations]')
[i m]
end

%------------ Stack it 20x --------------------------------------
% They stack the system to make it high dimensional.

Phi20 = Phi;
for i = 2:20
   Phi20 = cartProd(Phi20,Phi); 
end

% Construct linear tranform to untangle states, inputs, and updated states
order = [];
for i = 1:9
    order = [order i:9:180];
end
T = M_i(order,180);

% Untangle states, inputs, and updated states
Phi20 = T*Phi20;

%% ----------- Reachability Analysis using Phi20 --------------------

% Initial set represented as a hybrid zonotope
R{1} = VandM_2_HybZono(repmat(Apt{1},3,1),eye(10));
% All possible input
U = hybZono(.5*ones(60,1),zeros(60,1),diag(.5*ones(60,1)));
N = 10;

ZeroEye = [zeros(60,120)    eye(60)];
EyeZero = [eye(120)         zeros(120,60)];
tic
for i = 1:N
   startIt = toc;
   % Successor Set Identity
   R{i+1} = ZeroEye*andR(Phi20,cartProd(R{i},U),EyeZero);
   % Keeping track of computation times
   N_StepTime_TotalTime(i,:) = [i toc-startIt, toc];
end

% Display computation times
N_StepTime_TotalTime
% This takes roughly 1/10 of their reported times, but it seems like they
% do it 10x and we only do it once, so roughly the same time, assuming I am
% reading their example correctly.

%% ------------------SUBFUNCTIONS-------------------------

function [e] = e_i(i,n)
    e = zeros(1,n);
    e(i) = 1;
end

function [M] = M_i(inds,n)
    M = [];
    for i = 1:length(inds)
        e = e_i(inds(i),n);
        M = [M;e];
    end
end

function out = BooleanFunction(in,u)

    out = [ or(u(1),~xor(in(2),in(1)));
            ~xor(in(2),and(in(1),u(2)));
            ~and(in(3),~xor(u(2),u(3)))];

end

function Z = andR(obj,Y,R)
	Z = obj;
	
	if nargin < 3
		R = eye(size(Z.c,1));
	end
	
	% deal with Z not having an binary factors..
	if isempty(obj.Gb)
		obj.Gb = zeros(size(obj.c,1),0);
	end
	
	% use obj as a basis for the intersection
	Z.Gc = [ obj.Gc , zeros(size(obj.Gc,1),size(Y.Gc,2)) ];
	Z.Gb = [ obj.Gb , zeros(size(obj.Gb,1),size(Y.Gb,2)) ];
	Z.c = obj.c;
	% what dimensions are our constraints?
	ncz = max(size(Z.Ac,1),size(Z.Ab,1));
	ncy = max(size(Y.Ac,1),size(Y.Ab,1));
	% first combine the constraints
	Ac = [ obj.Ac , zeros(ncz,size(Y.Gc,2)) ;...
		zeros(ncy,size(obj.Gc,2)) , Y.Ac ];
	Ab = [ obj.Ab , zeros(ncz,size(Y.Gb,2)) ;...
			zeros(ncy,size(obj.Gb,2)) , Y.Ab ];
	% and intersection
	AcI = [ R*obj.Gc , -Y.Gc ];
	AbI = [ R*obj.Gb , -Y.Gb ];
	% put together
	Z.Ac = [ Ac ; AcI ];
	Z.Ab = [ Ab ; AbI ];
% 	Z.Ac = [ obj.Ac , zeros(ncz,size(Y.Gc,2)) ;...
% 		zeros(ncy,size(obj.Gc,2)) , Y.Ac ;...
% 		R*obj.Gc , -Y.Gc];
% 	if ~isempty(obj.Gb)
% 		Z.Ab = [ obj.Ab , zeros(ncz,size(Y.Gb,2)) ;...
% 			zeros(ncy,size(obj.Gb,2)) , Y.Ab ;...
% 			R*obj.Gb , -Y.Gb];
% 	else
% 		Z.Ab = [ obj.Ab , zeros(ncz,size(Y.Gb,2)) ;...
% 			zeros(ncy,size(obj.Gb,2)) , Y.Ab ];
% 	end
	Z.b = [ obj.b ; Y.b ; Y.c - R*obj.c ];

end		% end intersection