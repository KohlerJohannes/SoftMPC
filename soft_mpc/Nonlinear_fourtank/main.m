%This file runs the offline design and closed-loop simulations for the
% nonlinear four-tank system
%%
clear all
close all
clc
%% Offline design
Tank_offline;
load Results/w %load disturbance sequence
% initial conditions
t_0 = 0.0;   
%initial condition
x_init=0.45*ones(param.n,1);%[0.45;0.45;0.4;0.4];
u_init=2*ones(param.m,1);     
%simulation time
Tsim=param.t_switch(end);
boolean_compute=true; %runs solvers multiple times to make reliable estimates on comp. times; set to 1 
if boolean_compute
    compute_it=10;
else
    compute_it=1;
end
%% Define nonlinear MPC
%all MPC schems are defined next to each other, with objective obj, cost,
%constraints, and decision variables y; sing Casadi
import casadi.*
%% 
opts=struct;
opts.print_time=0;
opts.ipopt.print_level=0;
index_init=1:param.n;%index for initial condition
index_output_target=param.N*param.m+(param.N+1)*param.n+param.n+param.m+1:param.N*param.m+(param.N+1)*param.n+param.n+param.m+param.n/2;%index for desired output set in MPC
%proposed
y_relax=MX.sym('y_nom',(param.N*param.m+(param.N+1)*param.n)+param.n+param.m+param.n/2+param.n);
y_opt_relax=[repmat(x_init,param.N+1,1);repmat(u_init,param.N,1);x_init;u_init;param.y_des;x_init];
[solver_relax,con_lb_relax,con_ub_relax,lb_relax,ub_relax]=create_MPC_relax(param,y_relax,opts);
%proposed-robust
y_relax_robust=MX.sym('y_nom',(param.N*param.m+(param.N+1)*param.n)+param.n+param.m+param.n/2+param.n+1);
y_opt_relax_robust=[repmat(x_init,param.N+1,1);repmat(u_init,param.N,1);x_init;u_init;param.y_des;x_init;0];
[solver_relax_robust,con_lb_relax_robust,con_ub_relax_robust,lb_relax_robust,ub_relax_robust]=create_MPC_relax_robust(param,y_relax_robust,opts);
%nominal
y_nom=MX.sym('y_nom',(param.N*param.m+(param.N+1)*param.n)+param.n+param.m+param.n/2);
y_opt_nom=[repmat(x_init,param.N+1,1);repmat(u_init,param.N,1);x_init;u_init;param.y_des];
[solver_nom,con_lb_nom,con_ub_nom,lb_nom,ub_nom]=create_MPC_nom(param,y_nom,opts);
%soft
y_soft=MX.sym('y_soft',param.N*param.m+(param.N+1)*param.n+param.n+param.m+param.n/2+param.N*param.n);
y_opt_soft=[repmat(x_init,param.N+1,1);repmat(u_init,param.N,1);x_init;u_init;param.y_des;zeros(param.n*param.N,1)];
[solver_soft,con_lb_soft,con_ub_soft,lb_soft,ub_soft]=create_MPC_soft(param,y_soft,opts);
%initialize
x_soft(:,1)=x_init;x_nom(:,1)=x_init;x_relax(:,1)=x_init;x_relax_robust(:,1)=x_init;
c_nom_feasible=true;c_soft_feasible=true;
k_nom_infeasible=inf;k_soft_infeasible=inf;
t_relax=zeros(Tsim/param.delta,1);t_relax_robust=zeros(Tsim/param.delta,1);
t_soft=zeros(Tsim/param.delta,1);t_nom=zeros(Tsim/param.delta,1);
y_des_close=zeros(2,Tsim/param.delta);%initialize reference
tolerance=1e-6;%for checking feasiblity
%% Simulate closed-loop  
for k=1:Tsim/param.delta
%set reference 
y_des_close(:,k)=param.y_des(:,find(param.t_switch>=k*param.delta,1));
%% proposed
%candidate solution/warmstart
y_init_relax= [x_relax(:,k);y_opt_relax(param.n*2+1:param.n*(param.N+1));y_opt_relax(param.n*param.N+1:param.n*(param.N+1));                        y_opt_relax(param.n*(param.N+1)+param.m+1:param.n*(param.N+1)+param.m*param.N);...
    y_opt_relax(param.n*(param.N+1)+param.m*(param.N-1)+1:param.n*(param.N+1)+param.m*param.N);        y_opt_relax(end-(param.n*3/2+param.m)+1:end);            y_opt_relax(param.n+1:param.n*2)]; 
%set initial state
lb_relax(index_init)=x_relax(:,k);
ub_relax(index_init)=x_relax(:,k);
%set reference 
lb_relax(index_output_target)=y_des_close(:,k);
ub_relax(index_output_target)=y_des_close(:,k);
% Solve the NLP
t=tic;
for i=1:compute_it
res_relax = solver_relax('x0' , y_init_relax,... % solution guess
             'lbx', lb_relax,...           % lower bound on x
             'ubx', ub_relax,...           % upper bound on x
             'lbg', con_lb_relax,...           % lower bound on g
             'ubg', con_ub_relax);             % upper bound on g
end
t_relax(k)=toc(t)/compute_it;
%save result
y_opt_relax=full(res_relax.x); 
%simulate one step
u_relax_robust(:,k)=y_opt_relax_robust(param.n*(param.N+1)+1:param.n*(param.N+1)+param.m);
u_relax(:,k)=y_opt_relax(param.n*(param.N+1)+1:param.n*(param.N+1)+param.m);
r_relax(:,k)=y_opt_relax(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
%measure constraint violation 
con_violation_relax(k)=norm(min(x_relax(:,k)-param.x_min,0))+norm(max(x_relax(:,k)-param.x_max,0));
%check for infeasibility 
[~,ceq_relax_test]=nonlinearconstraints_relax(param,y_opt_relax);
if norm(ceq_relax_test)+sum(max(y_opt_relax-ub_relax,0))+sum(max(lb_relax-y_opt_relax,0)) >=tolerance
   error('relaxed infeasible') 
end
%% proposed-robust
%this is an improved formulation, which is covered in the extended online
%version (arxiv), but not included in the final manuscript
%candidate solution/warmstart
y_init_relax_robust= [x_relax_robust(:,k);y_opt_relax_robust(2*param.n+1:param.n*(param.N+1));y_opt_relax_robust(param.n*param.N+1:param.n*(param.N+1));  y_opt_relax_robust(param.n*(param.N+1)+param.m+1:param.n*(param.N+1)+param.m*param.N); ...
    y_opt_relax_robust(param.n*(param.N+1)+param.m*(param.N-1)+1:param.n*(param.N+1)+param.m*param.N); y_opt_relax_robust(end-1-(param.n*3/2+param.m)+1:end-1); y_opt_relax_robust(param.n+1:param.n*2);0]; 
%set initial state
lb_relax_robust(index_init)=x_relax_robust(:,k);
ub_relax_robust(index_init)=x_relax_robust(:,k);
%set reference 
lb_relax_robust(index_output_target)=y_des_close(:,k);
ub_relax_robust(index_output_target)=y_des_close(:,k);
% Solve the NLP
t=tic;
for i=1:compute_it
res_relax_robust = solver_relax_robust('x0' , y_init_relax_robust,... % solution guess
             'lbx', lb_relax_robust,...           % lower bound on x
             'ubx', ub_relax_robust,...           % upper bound on x
             'lbg', con_lb_relax_robust,...           % lower bound on g
             'ubg', con_ub_relax_robust);             % upper bound on g
end
t_relax_robust(k)=toc(t)/compute_it;
%save result
y_opt_relax_robust=full(res_relax_robust.x); 
%simulate one step
x_relax(:,k+1)=dynamic(x_relax(:,k),u_relax(:,k),param)+w(:,k);
x_relax_robust(:,k+1)=dynamic(x_relax_robust(:,k),u_relax_robust(:,k),param)+w(:,k);
r_relax_robust(:,k)=y_opt_relax_robust(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
%measure constraint violation 
con_violation_relax_robust(k)=norm(min(x_relax_robust(:,k)-param.x_min,0))+norm(max(x_relax_robust(:,k)-param.x_max,0));
%check for infeasibility 
[~,ceq_relax_robust_test]=nonlinearconstraints_relax_robust(param,y_opt_relax_robust);
if norm(ceq_relax_robust_test)+sum(max(y_opt_relax_robust-ub_relax_robust,0))+sum(max(lb_relax_robust-y_opt_relax_robust,0)) >=tolerance
   error('robust relaxed infeasible') 
end
%% Soft-MPC
if c_soft_feasible %first check if already became infeasible
    y_init_soft=[y_opt_soft(param.n+1:param.n*(param.N+1));y_opt_soft(param.n*param.N+1:param.n*(param.N+1));y_opt_soft(param.n*(param.N+1)+param.m+1:param.n*(param.N+1)+param.m*param.N);y_opt_soft(param.n*(param.N+1)+param.m*(param.N-1)+1:param.n*(param.N+1)+param.m*param.N);y_opt_soft(param.n*(param.N+1)+param.m*param.N+1:param.n*(param.N+1)+param.m*param.N+param.n*3/2+param.m);zeros(param.n*param.N,1)]; 
    lb_soft(index_init)=x_soft(:,k);
    ub_soft(index_init)=x_soft(:,k);
    lb_soft(index_output_target)=y_des_close(:,k);
    ub_soft(index_output_target)=y_des_close(:,k);
    t=tic;
    for i=1:compute_it
        res_soft = solver_soft('x0' , y_init_soft,... % solution guess
                     'lbx', lb_soft,...           % lower bound on x
                     'ubx', ub_soft,...           % upper bound on x
                     'lbg', con_lb_soft,...           % lower bound on g
                     'ubg', con_ub_soft);             % upper bound on g
    end
    t_soft(k)=toc(t)/compute_it;
    y_opt_soft=full(res_soft.x); 
    [c_soft_test,ceq_soft_test]=nonlinearconstraints_soft(param,y_opt_soft);%for soft constraint, need to check inequality constraints!
        if norm(ceq_soft_test)+sum(max(c_soft_test-con_ub_soft(size(c_soft_test,1)),0))+sum(max(y_opt_soft-ub_soft,0))+sum(max(lb_soft-y_opt_soft,0)) >=tolerance
           c_soft_feasible=false; 
           k_soft_infeasible=k;
        end
    u_soft(:,k)=y_opt_soft(param.n*(param.N+1)+1:param.n*(param.N+1)+param.m);
    x_soft(:,k+1)=dynamic(x_soft(:,k),u_soft(:,k),param)+w(:,k);%y_opt(param.n+1:2*param.n);
    r_soft(:,k)=y_opt_soft(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
    con_violation_soft(k)=norm(min(x_soft(:,k)-param.x_min,0))+norm(max(x_soft(:,k)-param.x_max,0));
end
%% nominal
if c_nom_feasible
    y_init_nom= [y_opt_nom(param.n+1:param.n*(param.N+1));y_opt_nom(param.n*param.N+1:param.n*(param.N+1));   y_opt_nom(param.n*(param.N+1)+param.m+1:param.n*(param.N+1)+param.m*param.N); y_opt_nom(param.n*(param.N+1)+param.m*(param.N-1)+1:param.n*(param.N+1)+param.m*param.N);y_opt_nom(end-(param.n*3/2+param.m)+1:end)]; 
    lb_nom(index_init)=x_nom(:,k);
    ub_nom(index_init)=x_nom(:,k);
    lb_nom(index_output_target)=y_des_close(:,k);
    ub_nom(index_output_target)=y_des_close(:,k);
    t=tic;
    for i=1:compute_it
        res_nom = solver_nom('x0' , y_init_nom,... % solution guess
                     'lbx', lb_nom,...           % lower bound on x
                     'ubx', ub_nom,...           % upper bound on x
                     'lbg', con_lb_nom,...           % lower bound on g
                     'ubg', con_ub_nom);             % upper bound on g
    end
   t_nom(k)=toc(t)/compute_it;
    y_opt_nom=full(res_nom.x); 
    [c_nom_test,ceq_nom_test]=nonlinearconstraints_nom(param,y_opt_nom);
    %for nominal, initial state constraint needs to be checked sperately
        if norm(ceq_nom_test)+sum(max(y_opt_nom-ub_nom,0))+sum(max(lb_nom-y_opt_nom,0))+sum(max(x_nom(:,k)-param.x_max,0))>=tolerance
           c_nom_feasible=false; 
           k_nom_infeasible=k;
        end
    u_nom(:,k)=y_opt_nom(param.n*(param.N+1)+1:param.n*(param.N+1)+param.m);
    x_nom(:,k+1)=dynamic(x_nom(:,k),u_nom(:,k),param)+w(:,k);%y_opt(param.n+1:2*param.n);
    r_nom(:,k)=y_opt_nom(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
    con_violation_nom(k)=norm(min(x_nom(:,k)-param.x_min,0))+norm(max(x_nom(:,k)-param.x_max,0));
end
end    
save('Results/Closed_loop')

%%
function cost = costfunction_nom(param,y)
    % Formulate the cost function to be minimized   
    cost = 0;   
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
% Build the cost by summing up the stage cost and the terminal cost
    for k=1:param.N
        x_k=x(param.n*(k-1)+1:param.n*k);
        u_k=u(param.m*(k-1)+1:param.m*k);
        cost = cost + runningcosts(x_k, u_k, r(1:param.n),r(param.n+1:param.n+param.m), param.Q, param.R);
    end
    %terminal cost:none
    offset=param.C*r(1:param.n)-y_target;
    cost=cost+param.alpha*offset'*offset; 
end
function cost = costfunction_soft(param,y)
    % Formulate the cost function to be minimized   
    cost = 0;   
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
y_slack=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n*param.N);
% Build the cost by summing up the stage cost and the terminal cost
    for k=1:param.N
        x_k=x(param.n*(k-1)+1:param.n*k);
        u_k=u(param.m*(k-1)+1:param.m*k);
        xi_k=y_slack(param.n*(k-1)+1:param.n*k);
        cost = cost + runningcosts(x_k, u_k, r(1:param.n),r(param.n+1:param.n+param.m), param.Q, param.R);
        %soft penalty
        cost = cost + xi_k'*param.Q_xi*xi_k;
    end
    %terminal cost - no (TEC) 
    %offset cost
    offset=param.C*r(1:param.n)-y_target; 
    cost=cost+param.alpha*offset'*offset;
end


function cost = costfunction_relax(param,y)
    % Formulate the cost function to be minimized   
    cost = 0;   
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
x_init=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n);
%last variable is nominal initial state
%-> make x trajectory nominal
x_nom=[x_init;x(param.n+1:end)];
x_real=x(1:param.n);
%cost initial state
 cost=cost+param.lambda*runningcosts(x_init, zeros(param.m,1),x_real ,zeros(param.m,1), param.P, param.R);
% Build the cost by summing up the stage cost and the terminal cost
    for k=1:param.N 
        x_k=x_nom(param.n*(k-1)+1:param.n*k);
        u_k=u(param.m*(k-1)+1:param.m*k);
        cost = cost + runningcosts(x_k, u_k, r(1:param.n),r(param.n+1:param.n+param.m), param.Q, param.R);
    end
    %terminal cost - no (TEC) 
    %offset cost
    offset=param.C*r(1:param.n)-y_target; 
    cost=cost+param.alpha*offset'*offset; 
end

function cost = costfunction_relax_robust(param,y)
    % Formulate the cost function to be minimized   
    cost = 0;   
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
x_init=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n);
slack=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n+1);
%last variable is nominal initial state
%-> make x trajectory nominal
x_nom=[x_init;x(param.n+1:end)];
%cost initial state
 cost=cost+param.lambda*slack^2;
% Build the cost by summing up the stage cost and the terminal cost
    for k=1:param.N 
        x_k=x_nom(param.n*(k-1)+1:param.n*k);
        u_k=u(param.m*(k-1)+1:param.m*k);
        cost = cost + runningcosts(x_k, u_k, r(1:param.n),r(param.n+1:param.n+param.m), param.Q, param.R);
    end
    %terminal cost:none
    %offset cost
    offset=param.C*r(1:param.n)-y_target;
    cost=cost+param.alpha*offset'*offset; 
end
 
function cost = runningcosts(x, u, x_eq, u_eq, Q, R)
    % Provide the running cost   
    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
    
end

%%
function [c, ceq] = nonlinearconstraints_soft(param,y) 
%    % Introduce the nonlinear constraints also for the terminal state  
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
y_slack=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n*param.N);
    c = [];
   ceq = [];
   % constraints along prediction horizon
    for k=1:param.N
        x_k=x((k-1)*param.n+1:k*param.n);
        xi_k=y_slack(param.n*(k-1)+1:param.n*k);
        x_new=x(k*param.n+1:(k+1)*param.n);        
        u_k=u((k-1)*param.m+1:k*param.m);
        %dynamic constraint
        ceqnew=x_new - dynamic(x_k, u_k,param);
        ceq = [ceq; ceqnew];
        %nonlinear constraints on state and input could be included here
        cnew=x_k+xi_k-param.x_max;
        c=[c;cnew];
        cnew=param.x_min-(x_k+xi_k);
        c=[c;cnew];
    end  
     %terminal set constraint - TEC
     %1. TEC
     ceqnew=x_new-r(1:param.n);
     ceq=[ceq;ceqnew];
     %2. setpoints
     ceqnew=fun_c(r(1:param.n),r(param.n+1:end),param);
     ceq=[ceq;ceqnew];
end
 

function [c, ceq] = nonlinearconstraints_nom(param,y) 
%    % Introduce the nonlinear constraints also for the terminal state  
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
    c = [];
   ceq = [];
   % constraints along prediction horizon
    for k=1:param.N
        x_k=x((k-1)*param.n+1:k*param.n);
        x_new=x(k*param.n+1:(k+1)*param.n);        
        u_k=u((k-1)*param.m+1:k*param.m);
        %dynamic constraint
        ceqnew=x_new - dynamic(x_k, u_k,param);
        ceq = [ceq; ceqnew]; 
    end  
     %terminal set constraint - TEC
     %1. TEC
     ceqnew=x_new-r(1:param.n);
     ceq=[ceq;ceqnew];
     %2. setpoints
     ceqnew=fun_c(r(1:param.n),r(param.n+1:end),param);
     ceq=[ceq;ceqnew];
end


function [c, ceq] = nonlinearconstraints_relax(param,y) 
%    % Introduce the nonlinear constraints also for the terminal state  
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
x_init=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n);
%last variable is nominal initial state
%-> make x trajectory nominal
x_nom=[x_init;x(param.n+1:end)];
x_real=x(1:param.n);
c = [];
   ceq = [];
   % constraints along prediction horizon
    for k=1:param.N
        x_k=x_nom((k-1)*param.n+1:k*param.n);
        x_new=x_nom(k*param.n+1:(k+1)*param.n);        
        u_k=u((k-1)*param.m+1:k*param.m);
        %dynamic constraint
        ceqnew=x_new - dynamic(x_k, u_k,param);
        ceq = [ceq; ceqnew];
        %nonlinear constraints on state and input could be included here
        %c=[c cnew];
    end  
     %terminal set constraint - TEC
     %1. TEC
     ceqnew=x_new-r(1:param.n);
     ceq=[ceq;ceqnew];
     %2. setpoints
     ceqnew=fun_c(r(1:param.n),r(param.n+1:end),param);
     ceq=[ceq;ceqnew];
end



function [c, ceq] = nonlinearconstraints_relax_robust(param,y) 
% Introduce the nonlinear constraints also for the terminal state  
x=y(1:param.n*(param.N+1));
u=y(param.n*(param.N+1)+1:param.n*(param.N+1)+param.N*param.m);
r=y(param.n*(param.N+1)+param.N*param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m);
y_target=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2);
x_init=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+1:param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n);
slack=y(param.n*(param.N+1)+param.N*param.m+param.n+param.m+param.n/2+param.n+1);
%last variable is nominal initial state
%-> make x trajectory nominal
x_nom=[x_init;x(param.n+1:end)];
x_real=x(1:param.n);
c = [];
   ceq = [];
   % constraints along prediction horizon
    for k=1:param.N
        x_k=x_nom((k-1)*param.n+1:k*param.n);
        x_new=x_nom(k*param.n+1:(k+1)*param.n);        
        u_k=u((k-1)*param.m+1:k*param.m);
        %dynamic constraint
        ceqnew=x_new - dynamic(x_k, u_k,param);
        ceq = [ceq; ceqnew];
    end  
     %terminal set constraint
     ceqnew=x_new-r(1:param.n);
     ceq=[ceq;ceqnew];
     %2. setpoints
     ceqnew=fun_c(r(1:param.n),r(param.n+1:end),param);
     ceq=[ceq;ceqnew];
     %initial state constraint
     cnew=sqrt(runningcosts(x_init, zeros(param.m,1),x_real,zeros(param.m,1), param.P, param.R)+1e-8)-(param.delta_relax+slack);%the factor 1e-8 provides a smooth approximation of the norm;
     %alternatively, this can be implemented as a quadratic constraints, to
     %avoid the 1e-8 slack for differentiability;
     %or a solver needs to be used that can handle conic constraitns
     c=[c;cnew];
end
 



%% create MPC
function [solver_relax,con_lb_relax,con_ub_relax,lb_relax,ub_relax]=create_MPC_relax(param,y_relax,opts)
import casadi.*
%obj_relax = MX(0);
obj_relax=costfunction_relax(param,y_relax);
[c_relax, ceq_relax] = nonlinearconstraints_relax(param,y_relax);
con_relax=[c_relax;ceq_relax];
con_bound_relax=zeros((param.N+2)*param.n,1);
con_lb_relax=[con_bound_relax];
con_ub_relax=[con_bound_relax];
lb_relax=[repmat(param.x_min,param.N+1,1);repmat(param.u_min,param.N,1);param.x_min;param.u_min;-inf(param.n/2,1);param.x_min];
ub_relax=[repmat(param.x_max,param.N+1,1);repmat(param.u_max,param.N,1);param.x_max;param.u_max;inf(param.n/2,1);param.x_max];
nlp_relax = struct('x', y_relax, 'f', obj_relax, 'g', con_relax);
solver_relax = nlpsol('solver', 'ipopt', nlp_relax,opts);  
end
function [solver_relax_robust,con_lb_relax_robust,con_ub_relax_robust,lb_relax_robust,ub_relax_robust]=create_MPC_relax_robust(param,y_relax_robust,opts);
import casadi.*
%obj_relax_robust = MX(0);
obj_relax_robust=costfunction_relax_robust(param,y_relax_robust);
[c_relax_robust, ceq_relax_robust] = nonlinearconstraints_relax_robust(param,y_relax_robust);
con_relax_robust=[c_relax_robust;ceq_relax_robust];
con_bound_relax_robust=zeros((param.N+2)*param.n,1);
con_lb_relax_robust=[-inf(1);con_bound_relax_robust];
con_ub_relax_robust=[0;con_bound_relax_robust];
lb_relax_robust=[repmat(param.x_min_robust,param.N+1,1);repmat(param.u_min,param.N,1);param.x_min_robust;param.u_min;-inf(param.n/2,1);param.x_min_robust;0];
ub_relax_robust=[repmat(param.x_max_robust,param.N+1,1);repmat(param.u_max,param.N,1);param.x_max_robust;param.u_max;inf(param.n/2,1);param.x_max_robust;inf];
nlp_relax_robust = struct('x', y_relax_robust, 'f', obj_relax_robust, 'g', con_relax_robust);
solver_relax_robust = nlpsol('solver', 'ipopt', nlp_relax_robust,opts);  
end
function [solver_nom,con_lb_nom,con_ub_nom,lb_nom,ub_nom]=create_MPC_nom(param,y_nom,opts);
import casadi.*
obj_nom=costfunction_nom(param,y_nom); 
[c_nom, ceq_nom] = nonlinearconstraints_nom(param,y_nom);
con_nom=[c_nom;ceq_nom];
con_bound_nom=zeros((param.N+2)*param.n,1);
con_lb_nom=[con_bound_nom];
con_ub_nom=[con_bound_nom];
lb_nom=[repmat(param.x_min,param.N+1,1);repmat(param.u_min,param.N,1);param.x_min;param.u_min;-inf(param.n/2,1)];
ub_nom=[repmat(param.x_max,param.N+1,1);repmat(param.u_max,param.N,1);param.x_max;param.u_max;inf(param.n/2,1)];
nlp_nom = struct('x', y_nom, 'f', obj_nom, 'g', con_nom);
solver_nom = nlpsol('solver', 'ipopt', nlp_nom,opts); 
end
function [solver_soft,con_lb_soft,con_ub_soft,lb_soft,ub_soft]=create_MPC_soft(param,y_soft,opts);
import casadi.*
obj_soft=costfunction_soft(param,y_soft);
[c_soft, ceq_soft] = nonlinearconstraints_soft(param,y_soft);
con_soft=[ceq_soft;c_soft];
con_bound_soft=zeros((param.N+2)*param.n,1);
con_lb_soft=[con_bound_soft;-inf(2*param.n*param.N,1)];
con_ub_soft=[con_bound_soft;zeros(2*param.n*param.N,1)];
lb_soft=[-inf(param.n*(param.N+1),1);repmat(param.u_min,param.N,1);param.x_min;param.u_min;-inf(param.n/2,1);-inf(param.n*param.N,1)];
ub_soft=[inf(param.n*(param.N+1),1); repmat(param.u_max,param.N,1); param.x_max;param.u_max;  inf(param.n/2,1);inf(param.n*param.N,1)];
nlp_soft = struct('x', y_soft, 'f', obj_soft, 'g', con_soft);
solver_soft = nlpsol('solver', 'ipopt', nlp_soft,opts); 
end