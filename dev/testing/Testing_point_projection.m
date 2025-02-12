%%
close all, clc, clear all

Gc = eye(2);
Gb = [2;0];
c = [0;0];
Z = hybZono(Gc, Gb, c, [], [], []);
y = randn(2,1);

figure
hold on, grid on, grid minor
plot(Z, 'b', .5)
plot(y(1), y(2), 'g.', 'MarkerSize',18)

axis([-3 3 -3 3])
axis equal

result = point_projection(Z,y);
x = result.x(1:Z.n);
feas_factors = result.x(Z.n+1:end);

plot(x(1), x(2), 'r.', 'MarkerSize', 18)


%%
clear all, close all, clc

% rng(3)
M = randn(2,5);
CZ1 = conZono([randn(2,2) [0;0]], [1;-1], randn(1,3), 0);
CZ2 = conZono([randn(2,2) [0;0]], [10;0], randn(1,3), 0);
% CZ3 = conZono(M, [7;-8], [randn(1,3) 0 0], 0);
% Z = hybZono(CZ1);
Z = union(CZ1, CZ2);
% Z = union(CZ1, union(CZ2, CZ3));
% Z = sharpHybZono(Z);
y = randn(2,1);

figure
hold on, grid on, grid minor
plot(Z, 'b', .5)
plot(y(1), y(2), 'g.', 'MarkerSize',18)
axis equal

result = point_projection(Z,y);
x = result.x(1:Z.n);
feas_factors = result.x((Z.n+1):end);

plot(x(1), x(2), 'r.', 'MarkerSize', 18)

%%
close all, clc, clear all

Gc = eye(2);
Gb = 2*eye(2);
c = [randn(2,1)];
Z = hybZono(Gc, Gb, c, [], [], []);
y = randn(2,1);

figure
hold on, grid on, grid minor
plot(Z, 'b', .5)
plot(y(1), y(2), 'g.', 'MarkerSize',18)

% axis([-3 3 -3 3])
axis equal

result = point_projection(Z,y);
x = result.x(1:Z.n);
feas_factors = result.x(Z.n+1:end);

plot(x(1), x(2), 'r.', 'MarkerSize', 18)
