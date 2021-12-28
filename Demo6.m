%% System Formulation
clear
close all
clc

A = [3 1 1;
     0 -5 0;
     0 2 3];
 
B = [1;2;1];

%% Decompose the System into Stable and Anti-stable Subsystems
[M,J,BLKS] = bdschur(A); % Jordan Decomposition

A_unstable = J(1:2,1:2);
          
B_z = inv(M)*B;
B_unstable = B_z(1:2);

%% Anti-stable Subsystem - Calculating NCR's Extremal Trajectories
boundary_R = @(t)  (2*((-1)^1)*expm(-A_unstable*(t-0))+((-1)^2)*eye(2))*inv(A_unstable)*B_unstable ;

x = [];
for t=0:0.1:100
    x = [x;round(boundary_R(t),3)'];
end

x = [x;-x];

plot(x(:,1),x(:,2))
xlabel('z1')
ylabel('z2')
title('NCR of Anti-stable Subsystem')

%% Overall System's NCR
figure
for i = -10:1:10
    dim_x = size(x);
    z = [x,kron(ones(dim_x(1),1),i)];
    x_full = z*M';   % Transform z-coordinate to x-coordinate
    plot3(x_full(:,1),x_full(:,2),x_full(:,3))
    hold on
end
xlabel('x1')
ylabel('x2')
zlabel('x3')
title('NCR of the Overall System')

%% Expressing CCLF as a "Look up Table"
k = boundary(x(:,1),x(:,2));
x = x(k,:);  % Updating x

x_dim = size(x);
x_num = x_dim(1);


p = [];
v = [];
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
x1 = 0.6328;
x2 = -9.178;
x3 = 1.878;
x = [x1 ; x2 ; x3];
z = inv(M)*x;   % Transform from x-coordimate to z-coordinate
z_unstable = [z(1);z(2)];  % Extract the anti-stable states [z1 z2]

lamda = 0.7;
umin = -1;
umax = 1;

t_sim = [];
x_sim = [];

for i=0:601
    u = u_CLF(z_unstable(1),z_unstable(2),A_unstable,B_unstable,umin,umax,lamda,F);
    dxdt = @(t,x) A*x+B*u; 
    [t,x_ode45] = ode45(dxdt,[0:0.01: 10],x);  
    x_sim = [x_sim;x_ode45(1:10,:)];
    t_sim = [t_sim;t(1:10)+0.1*i];
    x = x_ode45(11,:)';
    z = inv(M)*x;    % Transform from x-coordimate to z-coordinate
    z_unstable = [z(1);z(2)];  % Extract the anti-stable states [z1 z2]
end
hold on
plot3(x_sim(:,1),x_sim(:,2),x_sim(:,3))
hold off