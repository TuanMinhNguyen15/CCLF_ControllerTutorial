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

%% Simulation Results
%%%% CCLF vs LQR at initial state x = [-0.543 ; 0.85] %%%%
%%%%%% CCLF %%%%%%
hold on
x1 = -0.543;
x2 = 0.85;
z = [x1 ; x2];

lamda = 0.3;
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

%%%%%% LQR %%%%%%
Q = 2*eye(2);
R = 1;

K = lqr(A,B,Q,R);

x1 = -0.543;
x2 = 0.85;
z = [x1 ; x2];

umin = -1;
umax = 1;

t_sim = [];
z_sim = [];

for i=0:601
    u = -K*z;
    u = max(min(umax, u), umin);
    dzdt = @(t,z) A*z+B*u; 
    [t,z_ode45] = ode45(dzdt,[0:0.01: 10],z);  
    z_sim = [z_sim;z_ode45(1:10,:)];
    t_sim = [t_sim;t(1:10)+0.1*i];
    z = [z_ode45(11,1) ; z_ode45(11,2)];
end

plot(z_sim(:,1),z_sim(:,2))
axis([-1.5,1.5,-1.5,1.5]);

%%%% CCLF vs LQR at initial state x = [0.5 ; -0.5] %%%%
%%%%%% CCLF %%%%%%
hold on
x1 = 0.5;
x2 = -0.5;
z = [x1 ; x2];

lamda = 0.3;
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

%%%%%% LQR %%%%%%
Q = 2*eye(2);
R = 1;

K = lqr(A,B,Q,R);

x1 = 0.5;
x2 = -0.5;
z = [x1 ; x2];

umin = -1;
umax = 1;

t_sim = [];
z_sim = [];

for i=0:601
    u = -K*z;
    u = max(min(umax, u), umin);
    dzdt = @(t,z) A*z+B*u; 
    [t,z_ode45] = ode45(dzdt,[0:0.01: 10],z);  
    z_sim = [z_sim;z_ode45(1:10,:)];
    t_sim = [t_sim;t(1:10)+0.1*i];
    z = [z_ode45(11,1) ; z_ode45(11,2)];
end

plot(z_sim(:,1),z_sim(:,2))
axis([-1.5,1.5,-1.5,1.5]);
title('CCLF-Controller versus LQR-Controller')
xlabel('x1')
ylabel('x2')