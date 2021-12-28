%% System Formulation
clear 
close all
clc

A = [0 -0.5;
     1  1.5];
 
B = [0 1;1 0];

%% Subsystem 1 - Calculating NCR's Extremal Trajectories
A = [0 -0.5;
     1  1.5];
B1 = [0;1];

boundary_R = @(t)  (2*((-1)^1)*expm(-A*(t-0))+((-1)^2)*eye(2))*inv(A)*B1 ;

x1 = [];
for t=0:0.1:20
    x1 = [x1;round(boundary_R(t),3)'];
end

x1 = unique(x1,'rows');
x1 = [x1;-x1];


%% Subsystem 2 - Calculating NCR's Extremal Trajectories
A = [0 -0.5;
     1  1.5];
B2 = [1;0];

boundary_R = @(t)  (2*((-1)^1)*expm(-A*(t-0))+((-1)^2)*eye(2))*inv(A)*B2 ;

x2 = [];
for t=0:0.1:20
    x2 = [x2;round(boundary_R(t),3)'];
end

x2 = unique(x2,'rows');
x2 = [x2;-x2];


%% Plot NCR's of Subsystems 1&2
plot(x1(:,1),x1(:,2))
xlabel('x1')
ylabel('x2')
title("NCR of Sub-system 1")

figure 

plot(x2(:,1),x2(:,2))
xlabel('x1')
ylabel('x2')
title("NCR of Sub-system 2")

%% Minkowski Sum of subsystems 1&2 => NCR of the overall system
figure
P1 = Polyhedron('V',x1);
P1.minVRep();
P2 = Polyhedron('V',x2);
P2.minVRep();
P3 = P1.plus(P2);
P3.minVRep();
P3.plot()
xlabel('x1')
ylabel('x2')
title("NCR of Overall System")

%% Overall System - Finding NCR's Extremal Trajectories
x = [];
for theta = 0:0.01:2*pi
    d = [1;0];   % Shooting vector
    d = [cos(theta) -sin(theta);sin(theta) cos(theta)]*d;  % Rotating the shooting vector by the angle theta
    alpha = P3.shoot(d);  % Find the coefficient alpha such that alpha*d = a point on the NCR's boundary
    x_boundary = d*alpha.alpha; % Determine the point on the NCR's boundary
    x = [x;x_boundary'];
end

k = boundary(x(:,1),x(:,2));
x = x(k,:);  % Updating x

x_dim = size(x);
x_num = x_dim(1);


p = [];  % Data points vector
v = [];  % CCLF values vector
for mul = 0:0.005:1
    if mul == 0
        p = [p;[0 0]];
        v = [v;0];
    else
        p = [p;mul*x];
        v = [v;kron(ones(x_num,1),mul)];
    end
    
end

F = scatteredInterpolant(p,v);  % Expressing CCLF as a "Look-up Table" F

%% Simulation Results
figure
%%%% System's NCR %%%%
plot(x(:,1),x(:,2))
hold on

%%%% Trajectoriy starting at the point x = [-1;2] %%%%
x1 = -1;
x2 = 2;
z = [x1 ; x2];

lamda = 0.4;  % Convergence-rate control
umin = -1;
umax = 1;

t_sim = [];
z_sim = [];

for i=0:601
    u = u_CLF(z(1),z(2),A,B,umin,umax,lamda,F);  
    dzdt = @(t,z) A*z+B*u; 
    [t,z_ode45] = ode45(dzdt,[0:0.01: 10],z);  
    z_sim = [z_sim;z_ode45(1:10,:)];
    t_sim = [t_sim;t(1:10)+0.1*i];
    z = [z_ode45(11,1) ; z_ode45(11,2)];
end

plot(z_sim(:,1),z_sim(:,2))

%%%% Trajectoriy starting at the point x = [3;-2] %%%%
x1 = 3;
x2 = -2;
z = [x1 ; x2];

lamda = 0.4;
umin = -1;
umax = 1;

t_sim = [];
z_sim = [];

for i=0:601
    u = u_CLF(z(1),z(2),A,B,umin,umax,lamda,F);
    dzdt = @(t,z) A*z+B*u; 
    [t,z_ode45] = ode45(dzdt,[0:0.01: 10],z);  
    z_sim = [z_sim;z_ode45(1:10,:)];
    t_sim = [t_sim;t(1:10)+0.1*i];
    z = [z_ode45(11,1) ; z_ode45(11,2)];
end

plot(z_sim(:,1),z_sim(:,2))

%%%% Trajectoriy starting at the point x = [1;1] %%%%
x0 = 1;
y0 = 1;
z = [x0 ; y0];

lamda = 0.4;
umin = -1;
umax = 1;

t_sim = [];
z_sim = [];

for i=0:100
    u = u_CLF(z(1),z(2),A,B,umin,umax,lamda,F);
    dzdt = @(t,z) A*z+B*u; 
    [t,z_ode45] = ode45(dzdt,[0:0.01: 10],z);  
    z_sim = [z_sim;z_ode45(1:10,:)];
    t_sim = [t_sim;t(1:10)+0.1*i];
    z = [z_ode45(11,1) ; z_ode45(11,2)];
end

plot(z_sim(:,1),z_sim(:,2))
axis([-4,4,-3,3]);
title('CCLF Controller')
xlabel('x1')
ylabel('x2')

