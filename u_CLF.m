function u = u_CLF(x0,y0,A,B,umin,umax,lamda,F)
    dim_B = size(B);
    n = dim_B(1);
    m = dim_B(2);
    dVdz = dVdzCal(x0,y0,F,0.000000001,15);
    Aineq = [dVdz*B -1];
    bineq = -lamda*F(x0,y0)-dVdz*A*[x0;y0];
    Aeq = [];
    beq = [];
    lb = [kron(ones(m,1),umin);0];
    ub = [kron(ones(m,1),umax);inf];
    H = diag([ones(1,m),0]);
    f = [zeros(1,m)  1000];
    sol = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub);
    u = sol(1:m);
end