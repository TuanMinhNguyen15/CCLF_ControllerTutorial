%% System Formulation
clear
close all
clc

A = [0 -0.5;
     1 1.5];
 
B = [0;
     -1];
 
%% Calculating NCR's Extremal Trajectories
boundary_R = @(t)  (2*((-1)^1)*expm(-A*(t-0))+((-1)^2)*eye(2))*inv(A)*B ;

x = [];
for t=0:0.1:100
    x = [x;round(boundary_R(t),3)'];
end

x = unique(x,'rows');
x = [x;-x];

plot(x(:,1),x(:,2))
xlabel('x1');
ylabel('x2');
title("NCR's Boundary")

%% Expressing CCLF as a "Look up Table"
x_dim = size(x);
x_num = x_dim(1);

p = [];  % Data points vector
v = [];  % CCLF values vector
for mul = 0:0.01:1 % 0:0.005:1
    if mul == 0
        p = [p;[0 0]];
        v = [v;0];
    else
        p = [p;mul*x];
        v = [v;kron(ones(x_num,1),mul)];
    end
    
end

F = scatteredInterpolant(p,v); % Expressing CCLF as a "Look-up Table" F

%% Plotting CCLF
figure

V = F(p);  % Evaluating CCLF values V at points p via the "Look-up Table" F

plot3(p(:,1),p(:,2),V)
title('CCLF Plot')
xlabel('x1')
ylabel('x2')
zlabel('V')