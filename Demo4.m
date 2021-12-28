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

x = unique(x,'rows');  % Filter out the duplicated boundary points
x = [x;-x];            % x = [Bunch-2 Trajectories ; Bunch-1 Trajectories]

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

%% Plotting CCLF V and its Euler Approximation of dV/dx at the point x = [0.6 ; -0.3]

%%%% Plotting CCLF V %%%%
V = F(p);  % Evaluating CCLF values V at points p via the "Look-up Table" F

plot3(p(:,1),p(:,2),V)
title('CCLF Plot')
xlabel('x1')
ylabel('x2')
zlabel('V')

%%%% Euler Approximation %%%%
d = 0.001; % Grid distance
n = 8;     % Number of grid points in each axis

x1 = .6;
x1max = x1+d;
x1min = x1-d;
x1_grid = x1min:2*d/n:x1max;  % Grid points in the x1-axis

x2 = -.3;
x2max = x2+d;
x2min = x2-d;
x2_grid = x2min:2*d/n:x2max;  % Grid points in the x2-axis

Pq = [];   % Grid points matrix

for i=1:length(x1_grid)
    for j = 1:length(x2_grid)
        Pq = [Pq;x1_grid(i) , x2_grid(j)];
    end
end


Vq = F(Pq);  % Evaluating CCLF values Vq at the grid points in matrix Pq
mdl = fitlm(Pq,Vq);  % Fit a linear model to the sample points [Pq ; Vq]
result_table = mdl.Coefficients;  % Extracting the coefficients of the linear model
result_array = table2array(result_table);  % Convert the table result into an array

dVdx = @(x1,x2) result_array(2,1)*x1 + result_array(3,1)*x2 + result_array(1,1);  % Construct a hyperplane from the linear model
hold on
fsurf(dVdx)
plot3(x1,x2,F(x1,x2),'-o','Color','b','MarkerSize',10,...   % Plot the hyperplane
    'MarkerFaceColor','#D9FFFF')
axis([-2 2 -2 2 0 1])
hold off