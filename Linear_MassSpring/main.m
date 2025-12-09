%%
% This file runs the offline design and closed-loop simulations for the
% linear mass-srping-damper example
% The code for the Soft-T implementation is taken from [Wabersich et al. ]
%%
clear all
close all
clc
%
n=2;m=1;%state/input dimension
%dynamic equations
m=1;k=1;d=0.1;
A_c=[0,1;-k/m,-d/m];
B_c=[0;1];
h=0.05;%sampling time
[A,B]=c2d(A_c,B_c,h);%discretization
%cost
Q=diag([1,0.1]);R=0.2;
%state constraints: Ax*x<=bx
Ax = [1 0; -1 0; 0, 1; 0, -1];
bx = [1; 1;1;1];
%input constraints: Au*u<= bu
Au = [1;-1];
bu = 2*[1; 1];
% system info
n=size(Ax,2);
m=size(Au,2);
nax=size(Ax,1);
nau=size(Au,1);
% define LQR controller
[K_LQR,P_LQR]=dlqr(A,B,Q,R);
K_LQR=-K_LQR;
%Lyapunov equaiont
P_OL=dlyap(A',eye(n));
rho_delta=max(eig((A'*P_OL*A)/P_OL));
%terminal set with & witout state constraints
[~,G_x,g_x]=maxInvSet(A+B*K_LQR,[Ax;Au*K_LQR],[bx;bu]);
[~,G_u,g_u]=maxInvSet(A+B*K_LQR,[Au*K_LQR],[bu]);
x_symb=sdpvar(2,1);
constraint_x=[G_x*x_symb <=g_x];
%O_Inf_x=YSet(x_symb(:),constraint_x);
constraint_u=[G_u*x_symb <=g_u];
constraint_u_bounded=[G_u;eye(n);-eye(n)]*x_symb <=[g_u;1e2*ones(2*n,1)];
%O_Inf_u=YSet(x_symb(:),constraint_u_bounded);
%% cost
ell = @(x,u) x' * Q * x + u' * R * u;
c_xi=1e5;
lambda=1e3;
Q_xi = c_xi*eye(nax);
ell_xi = @(x,u,xi) ell(x,u)+(xi)'*Q_xi*(xi);
ell_f_LQR =@(x) x'*P_LQR*x;
ell_f_LQR_xi =@(x,xi) ell_f_LQR(x)+(xi)'*Q_xi*(xi);
%terminal cost for Soft-G 
Q_upper=Q+Ax'*Q_xi*Ax;
P_global=dlyap(A',Q_upper);
ell_f_global=@(x) x'*P_global*x;
%proposed initial cost
V_delta=@(x) lambda*x'*P_OL*x;
N=10; %prediction horizon
%% offline design for Soft-T
%taken from (Wabersich et al.)
%compute dual variables mu
mu=zeros(nax,1);
for i=1:nax
    options = optimoptions('linprog','Display','none');
    [~,mu(i)]=linprog(g_u,[],[],G_u',Ax(i,:)',zeros(size(g_u)),[],[],options);
end
%% Set up MPCs
yalmip_optimizer_nominal=create_mpc_nominal(N, A, B, Ax, bx, Au, bu,G_x,g_x, ell,ell_f_LQR);
yalmip_optimizer_proposed=create_mpc_proposed(N, A, B, Ax, bx, Au, bu,G_x,g_x, ell,ell_f_LQR,V_delta);
yalmip_optimizer_soft_P=create_mpc_soft_terminalhard(N, A, B, Ax, bx, Au, bu,G_x,g_x, ell_xi,ell_f_LQR);
yalmip_optimizer_soft_T=create_mpc_soft_terminalRelax(N, A, B, Ax, bx, Au, bu,G_u,g_u, ell_xi,ell_f_LQR_xi,mu);
yalmip_optimizer_soft_G=create_mpc_global(N, A, B, Ax, bx, Au, bu, ell_xi,ell_f_global);
%% Closed-loop Simulations
Tsim=300;
for mod =[1,3,4]%cycle 
 clear x_traj_nominal u_traj_nominal obj_traj_nominal cost_traj_nominal constr_traj_nominal T_nominal... 
       x_traj_soft_P u_traj_soft_P obj_traj_soft_P cost_traj_soft_P constr_traj_soft_P T_soft_P...
       x_traj_soft_T u_traj_soft_T obj_traj_soft_T cost_traj_soft_T constr_traj_soft_T T_soft_T...
       x_traj_soft_G u_traj_soft_G obj_traj_soft_G cost_traj_soft_G constr_traj_soft_G T_soft_G...
       x_traj_prop u_traj_prop obj_traj_prop cost_traj_prop constr_traj_prop T_prop
x0=[0.832;1];%initial conditions picked by hand on boundary of feasible sets for different approaches
if mod==1
    x0=1*x0;
elseif mod==2
    x0=1.017*x0;
elseif mod==3
    x0=1.520*x0;
 elseif mod==4
    x0=4*x0;     
 else 
     error('undefned mod/initial condition')
 end
 %select which MPCs to run based on feasibility
if mod==1
[x_traj_nominal, u_traj_nominal, obj_traj_nominal, cost_traj_nominal,constr_traj_nominal] = simulation(yalmip_optimizer_nominal,x0, A, B,Tsim, ell,Ax,bx);
end
if mod==1|| mod==2
[x_traj_soft_P, u_traj_soft_P, obj_traj_soft_P, cost_traj_soft_P,constr_traj_soft_P] = simulation(yalmip_optimizer_soft_P,x0, A, B,Tsim, ell,Ax,bx);
end
if mod==1|| mod==2||mod==3
[x_traj_soft_T, u_traj_soft_T, obj_traj_soft_T, cost_traj_soft_T,constr_traj_soft_T] = simulation(yalmip_optimizer_soft_T,x0, A, B,Tsim, ell,Ax,bx);
end
[x_traj_soft_G, u_traj_soft_G, obj_traj_soft_G, cost_traj_soft_G,constr_traj_soft_G] = simulation(yalmip_optimizer_soft_G,x0, A, B,Tsim, ell,Ax,bx);
[x_traj_prop, u_traj_prop, obj_traj_prop, cost_traj_prop,constr_traj_prop] = simulation(yalmip_optimizer_proposed,x0, A, B,Tsim, ell,Ax,bx);
save(['results/traj_mod_' num2str(mod)])
end
%% Run timming experiments
%% Gridd state space
Delta_x=2e-2;
[X1,X2]=meshgrid(-1:Delta_x:1,-1:Delta_x:1);
t_nom=0*X1;t_prop=t_nom;t_G=t_nom;t_P=t_nom;t_T=t_nom;%initialize by 0
feasible=false(size(X1));
%time is measured using "timeit", which calls the function multiple times
%and takes median to avoid fluctuations from process and overhead
for i=1:size(X1,1)
    for j=1:size(X1,2)
       x_temp=[X1(i,j);X2(i,j)];
       [res,flag] = yalmip_optimizer_nominal(x_temp); % call mpc
        if(flag==0||flag==3||flag==5)
           feasible(i,j)=true;
           timer_nominal=@()yalmip_optimizer_nominal(x_temp);%make function
           t_nom(i,j)=timeit(timer_nominal);
           timer_proposed=@() yalmip_optimizer_proposed(x_temp);
           t_prop(i,j)=timeit(timer_proposed);
           timer_soft_G=@() yalmip_optimizer_soft_G(x_temp);
           t_G(i,j)=timeit(timer_soft_G);
           timer_soft_T=@() yalmip_optimizer_soft_T(x_temp);
           t_T(i,j)=timeit(timer_soft_T);
           timer_soft_P=@() yalmip_optimizer_soft_P(x_temp);
           t_P(i,j)=timeit(timer_soft_P);
        end
    end
end
%%
save('results/comp_time.mat','t_P','t_G','t_nom','t_prop','feasible')
%% Auxillary functions
function yalmip_optimizer = create_mpc_soft_terminalRelax(N, A, B, Ax, bx, Au, bu,G,g, ell_xi,ell_f_LQR_xi,mu)
%soft-constraint MPC; relaxed terminal; 
    % for better readability
    n = size(A,1);
    m = size(B,2);
    nax = size(Ax, 1);
    % define optimization variables
    U = sdpvar(m,N-1,'full');
    X = sdpvar(n,N,'full');
    Xi = sdpvar(nax,N,'full');
    x0 = sdpvar(n,1,'full');
    alpha=sdpvar(1,1,'full');
    objective = 0;
    constraints = [X(:,1) == x0]; % initial condition
    constraints = [constraints, Xi(:,N) >= 0]; % terminal slack variables greater zero
    for k = 1:N-1
        constraints = [constraints, Xi(:,k) >= 0]; % slack variables greater zero
        constraints = [constraints, X(:,k+1) == A * X(:,k) + B * U(:,k)]; % dynamics constraint
        constraints = [constraints, Au * U(:,k) <= bu]; % hard input constraint
        constraints = [constraints, Ax * X(:,k) <= bx + Xi(:,k) + Xi(:,N)]; % soft state constraint
        objective = objective + ell_xi(X(:,k),U(:,k),Xi(:,k)+Xi(:,N));
        %objective + ell(X(:,k),U(:,k)) + c1 * ones(1,nax) * Xi(:,k) + c1 * ones(1,nax) * Xi(:,N); % stage cost
    end
    % terminal objective
    objective = objective + ell_f_LQR_xi(X(:,N),Xi(:,N)); % terminal slack
     %terminal set constraint
     constraints = [constraints, G * X(:,N) <= g*alpha];
     constraints = [constraints, alpha >= 0];
     constraints = [constraints, alpha <= 1];
     constraints = [constraints, alpha*mu<=bx+Xi(:,N)];
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');
    yalmip_optimizer = optimizer(constraints,objective,ops,x0,{objective, U, X, Xi, alpha});
end
%%
function yalmip_optimizer = create_mpc_soft_terminalhard(N, A, B, Ax, bx, Au, bu,G,g, ell_xi,ell_f_LQR)
%soft-constraint MPC; hard terminal penalty based on G_x
    % for better readability
    n = size(A,1);
    m = size(B,2);
    nax = size(Ax, 1);
    % define optimization variables
    U = sdpvar(m,N-1,'full');
    X = sdpvar(n,N,'full');
    Xi = sdpvar(nax,N-1,'full');
    x0 = sdpvar(n,1,'full');
    objective = 0;
    constraints = [X(:,1) == x0]; % initial condition
    for k = 1:N-1
        constraints = [constraints, Xi(:,k) >= 0]; % slack variables greater zero
        constraints = [constraints, X(:,k+1) == A * X(:,k) + B * U(:,k)]; % dynamics constraint
        constraints = [constraints, Au * U(:,k) <= bu]; % hard input constraint
        constraints = [constraints, Ax * X(:,k) <= bx + Xi(:,k)]; % soft state constraint
        objective = objective + ell_xi(X(:,k),U(:,k),Xi(:,k));
        %objective + ell(X(:,k),U(:,k)) + c1 * ones(1,nax) * Xi(:,k) + c1 * ones(1,nax) * Xi(:,N); % stage cost
    end
    % terminal objective
    objective = objective + ell_f_LQR(X(:,N)); 
     %terminal set constraint
     constraints = [constraints, G * X(:,N) <= g];
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');
    yalmip_optimizer = optimizer(constraints,objective,ops,x0,{objective, U, X, Xi});
end
%%
function yalmip_optimizer = createOL(A, B)
%
%proposed MPC with relaxed initial penalty
    % for better readability
    n = size(A,1);
    m = size(B,2);
    % define optimization variables
    U = sdpvar(m,1,'full');
    x0 = sdpvar(n,1,'full');
    objective = U'*U;
    constraints = []; 
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');
    yalmip_optimizer = optimizer(constraints,objective,ops,x0,{objective, U});
end
%%
function yalmip_optimizer = create_mpc_nominal(N, A, B, Ax, bx, Au, bu,G,g, ell,ell_f_LQR)
%nominal MPC
    n = size(A,1);
    m = size(B,2);
    % define optimization variables
    U = sdpvar(m,N-1,'full');
    X = sdpvar(n,N,'full');
    x0 = sdpvar(n,1,'full');
    objective = 0;
    constraints = [X(:,1) == x0]; % initial condition
    for k = 1:N-1
        constraints = [constraints, X(:,k+1) == A * X(:,k) + B * U(:,k)]; % dynamics constraint
        constraints = [constraints, Au * U(:,k) <= bu]; % hard input constraint
        constraints = [constraints, Ax * X(:,k) <= bx]; % hard state constraint
        objective = objective + ell(X(:,k),U(:,k));
    end
    % terminal objective
    objective = objective + ell_f_LQR(X(:,N)); 
    constraints = [constraints, G * X(:,N) <= g];
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');
    yalmip_optimizer = optimizer(constraints,objective,ops,x0,{objective, U, X});
end
%%
function yalmip_optimizer = create_mpc_proposed(N, A, B, Ax, bx, Au, bu,G,g, ell,ell_f_LQR,V_delta)
%proposed MPC with relaxed initial penalty
    % for better readability
    n = size(A,1);
    m = size(B,2);
    nax = size(Ax, 1);
    % define optimization variables
    U = sdpvar(m,N-1,'full');
    X = sdpvar(n,N,'full');
    x0 = sdpvar(n,1,'full');
    objective = 0;
    %constraints = [X(:,1) == x0]; % initial condition
    objective=objective+V_delta(X(:,1)-x0);
    constraints=[];
    for k = 1:N-1
        constraints = [constraints, X(:,k+1) == A * X(:,k) + B * U(:,k)]; % dynamics constraint
        constraints = [constraints, Au * U(:,k) <= bu]; % hard input constraint
        constraints = [constraints, Ax * X(:,k) <= bx]; % hard state constraint
        objective = objective + ell(X(:,k),U(:,k));
    end
    % terminal objective
    objective = objective + ell_f_LQR(X(:,N)); 
    constraints = [constraints, G * X(:,N) <= g];
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');
    yalmip_optimizer = optimizer(constraints,objective,ops,x0,{objective, U, X});
end
%%
function yalmip_optimizer = create_mpc_global(N, A, B, Ax, bx, Au, bu, ell_xi,ell_f_global)
%soft-constraint MPC;no terminal set, globally valid terminal cst
    n = size(A,1);
    m = size(B,2);
    nax = size(Ax, 1);
    % define optimization variables
    U = sdpvar(m,N-1,'full');
    X = sdpvar(n,N,'full');
    Xi = sdpvar(nax,N-1,'full');
    x0 = sdpvar(n,1,'full');
    objective = 0;
    constraints = [X(:,1) == x0]; % initial condition
    for k = 1:N-1
        constraints = [constraints, Xi(:,k) >= 0]; % slack variables greater zero
        constraints = [constraints, X(:,k+1) == A * X(:,k) + B * U(:,k)]; % dynamics constraint
        constraints = [constraints, Au * U(:,k) <= bu]; % hard input constraint
        constraints = [constraints, Ax * X(:,k)-Xi(:,k) <= bx ]; % soft state constraint
        objective = objective + ell_xi(X(:,k),U(:,k),Xi(:,k));
    end
    % terminal objective
    objective = objective + ell_f_global(X(:,N)); 
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');
    yalmip_optimizer = optimizer(constraints,objective,ops,x0,{objective, U, X, Xi});
end
%%

function [x_traj, u_traj, obj_traj, cost_traj,constr_traj] = simulation(yalmip_optimizer, x0, A, B, Tsim, ell,Ax,bx)
    u_traj = [];
    x_traj = [x0];
    obj_traj = [];
    cost_traj = 0; 
    
    for k = 1:1:Tsim
        [res,flag] = yalmip_optimizer(x_traj(:,end)); % call mpc
        assert(flag==0||flag==3||flag==5);%0:sucess, 3:max-iteration, 5=lack of progress, 
        obj_traj(k) = res{1};
        u_sol_mpc{k} = res{2};
        u_traj(:,end+1) = u_sol_mpc{k}(:,1);
        x_traj(:,end+1) = A * x_traj(:,end) + B * u_traj(:,end);
        %cost
        cost_traj(k) = ell(x_traj(:,k),u_traj(:,k));
        %constraints
        constr_traj(:,k)  = max(Ax*x_traj(:,k)-bx,0);
    end
    x_traj(:,end)=[];
end 

function [O_Inf,G,g]=maxInvSet(A,H,h)
% computes the maximal output admissible set for a discrete time linear
% system based on [Gilbert, Tan: Linear Systems with State and Control
% Constraints: The Theory and Application of Maximal Output Admissible Sets,
% IEEE Transactions on Automatic Control, vol.36, No. 9, 1991
% autonomous system: x+ = Ax
% constraints:       H*x-h <= 0
% MOAS O_inf defined by G*x <= g
options = optimset('Display','off');
m=length(h(:,1));
notFinished=1;
fmax=-inf;
h_new=h;
H_new=H;

while(notFinished)
 for i=1:m 
    [~,fval]=linprog(-H_new(end-m+i,:)*A,H_new,h_new,[],[],[],[],[],options);
    fmax=max(-fval-h(i),fmax);
 end
 if(fmax<=0)
     notFinished=0;
 else
     fmax=-inf;   
     H_new=[H_new; H_new(end-m+1:end,:)*A];
     h_new=[h_new;h_new(end-m+1:end)];   
 end
end
G=H_new;
g=h_new;
O_Inf = Polyhedron('A',G,'b',g);
%O_Inf.minHRep();
%O_Inf.minVRep();
end