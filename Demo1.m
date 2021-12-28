%% System Formulation
clear
close all
clc

A = [0 -0.5;
     1 1.5];
 
B = [0;
     -1];
 
%% Calculating NCR's Extremal Trajectories 
boundary_R = @(t)  (2*((-1)^1)*expm(-A*(t-0))+((-1)^2)*eye(2))*inv(A)*B ;  % Extremal trajectories formula

x = [];   % Initialize vector x which contains points on the NCR's boundary
for t=0:0.01:100
    x = [x;round(boundary_R(t),3)'];    % Iteratively adding boundary points to vector x
end

x = unique(x,'rows');  % Filter out the duplicated boundary points


plot(-x(:,1),-x(:,2))  % Bunch-1 Trajectories
hold on
plot(x(:,1),x(:,2))    % Bunch-2 Trajectories
hold off
xlabel('x1');
ylabel('x2');
title("NCR's Boundary")
legend('Bunch-1 Trajectories','Bunch-2 Trajectories');

