%offline computations
clear all 
close all
clc 
%% 
param=Tank_param();%define model
syms x [param.n,1]; 
syms u [param.m,1];
f_c_func = matlabFunction(fun_c(x,u,param),'Vars',{x,u});
% Define cont.-time derivatives
A = jacobian(fun_c(x,u,param),x);
B = jacobian(fun_c(x,u,param),u);
A_func = matlabFunction(A,'Vars',{x,u});
B_func = matlabFunction(B,'Vars',{x,u}); 
%analytical bounds for matrix P (see Appendix D)
p_3_min=param.c_3*sqrt(max(param.x_max_real)/min(param.x_min_real))/(4*param.c_1);
p_4_min=param.c_4*sqrt(max(param.x_max_real)/min(param.x_min_real))/(4*param.c_2);
P=param.P;
if P(3,3)<=p_3_min ||P(4,4)<=p_4_min
   error('Lyapunov does not hold (matrix P chosen incorrectly)') 
end
%% Compute contraction rate
rho_c=inf;
n_grid=2;%since functions are monotone, just vertices suffices; but a finer grid acn be chosen to double check
for x1=linspace(param.x_min_real(1),param.x_max_real(1),n_grid)
for x2=linspace(param.x_min(2),param.x_max(2),n_grid)
for x3=linspace(param.x_min(3),param.x_max(3),n_grid)
for x4=linspace(param.x_min(4),param.x_max(4),n_grid)
             A_temp=A_func([x1;x2;x3;x4],zeros(param.m,1));
             norm(A_temp+0.5*[param.c_1/sqrt(x1),0,-param.c_3/sqrt(x3),0;...
                         0,param.c_2/sqrt(x2),0,-param.c_4/sqrt(x4);...
                         0,0,param.c_3/sqrt(x3),0;...
                         0,0,0,param.c_4/sqrt(x4)]);                 
             rho_c=min(-max(eig((P*A_temp+A_temp'*P)/P))/2,rho_c);
             
end
end
end 
end 
% compute discrete-time contraction + max distu
param.P=P;
param.rho=exp(-rho_c*param.delta);%discrete-time contraction rate
param.w_bar=(1-param.rho)/sqrt(norm(P,2))*param.delta_relax; 