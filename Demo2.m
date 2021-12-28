%% System Formulation
clear
close all
clc

A = [0.2   1   0;
     0     0.2 0;
     0     0   0.4];
 
B = [1;1;1];

%% Calculating NCR's Extremal Trajectories - Bunch 1
Boundary_R = @(t,t1,t2) (2*((-1)^1)*expm(-A*(t-t1))+2*((-1)^2)*expm(-A*(t-t2))+((-1)^3)*eye(3))*inv(A)*B;  % Extremal trajectories formula

for t2 = 1:50
    x = [];
    t1 = 0;
    
    for t = t2:0.1:100
        x = [x;Boundary_R(t,t1,t2)'];
    end

    %plot3(x(:,1),x(:,2),x(:,3))
    hold on
    plot3(-x(:,1),-x(:,2),-x(:,3))
    
    title('Bunch-1 Trajectories');
    xlabel('x1');
    ylabel('x2');
    zlabel('x3');
end


%% Calculating NCR's Extremal Trajectories - Bunch 2
figure

for t2 = 1:50
    x = [];
    t1 = 0;
    
    for t = t2:0.1:100
        x = [x;Boundary_R(t,t1,t2)'];
    end

    plot3(x(:,1),x(:,2),x(:,3))
    hold on
    %plot3(-x(:,1),-x(:,2),-x(:,3))
    
    title('Bunch-2 Trajectories');
    xlabel('x1');
    ylabel('x2');
    zlabel('x3');
end
%% Calculating NCR's Extremal Trajectories - Overall
figure

for t2 = 1:50
    x = [];
    t1 = 0;
    
    for t = t2:0.1:100
        x = [x;Boundary_R(t,t1,t2)'];
    end

    plot3(x(:,1),x(:,2),x(:,3))
    hold on
    plot3(-x(:,1),-x(:,2),-x(:,3))
    
    title('Overall Extremal Trajectories');
    xlabel('x1');
    ylabel('x2');
    zlabel('x3');
end



