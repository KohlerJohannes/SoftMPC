
function param=Tank_param()
%% Example from [46, Limon2018, Sec VI] 
param.n = 4; 
param.m = 2;
%set parameters (SI units)
param.g=9.81;
param.S=0.06;
param.a_1=1.2938e-4;
param.a_2=1.5041e-4;
param.a_3=1.0208e-4;
param.a_4=9.3258e-5;
param.gamma_a=0.3;
param.gamma_b=0.4;
%also set constrains
param.u_min=[0;0];
param.u_max=[3.6;4];
%real state constraints
param.x_min=0.2*ones(4,1);
param.x_max=[1.36;1.36;1.30;1.30];
%maximal minimal state for offline
param.x_min_real=2e-2*ones(4,1);
param.x_max_real=2*ones(4,1);

param.V_max=0.2226;
%sampling time
param.delta=10;
%compute directly relevant constants
param.c_1=param.a_1/param.S*sqrt(2*param.g);
param.c_2=param.a_2/param.S*sqrt(2*param.g);
param.c_3=param.a_3/param.S*sqrt(2*param.g);
param.c_4=param.a_4/param.S*sqrt(2*param.g);
param.c_u1=param.gamma_a/(3600*param.S);
param.c_u2=param.gamma_b/(3600*param.S);
param.c_u3=(1-param.gamma_b)/(3600*param.S);
param.c_u4=(1-param.gamma_a)/(3600*param.S);

%equlibrium level:
param.x_0=[0.6702;0.6549;0.5435;0.5887];
param.u_0=[1.63;2];
%desired setpoints y=x_1,x_2 
param.y_des=[param.x_max(1);param.x_max(2)];
param.t_switch=7e3; 
%output map
param.C=[eye(2),zeros(2)];
%tec-> same constraints artifical
param.yr_min=param.x_min(1:2);%0.43;
param.yr_max=param.x_max(1:2);
%cost:
param.Q=eye(param.n);
param.R=1e-2*eye(param.m);
%weight offset cost
param.alpha=1e4;
%horizon
param.N=10; 
%disturbance
%param.Sigma=eye(param.n)*1e-3;
param.w_max=5e-2;
%param.disturbances_t=[5e2,2.5e3,4.5e3];
param.disturbances_t=[5e2,3.5e3,6.5e3];
%soft 
param.Q_xi=1e4*eye(param.n);
param.P=diag([1,1,2,2]);
param.lambda=1e5;
%set size of RPI set \bar{delta}
param.delta_relax=5e-2;
%tightened constraints
param.x_max_robust=param.x_max-param.delta_relax./sqrt(diag(param.P));
param.x_min_robust=param.x_min+param.delta_relax./sqrt(diag(param.P));
%[1500,3000];
end
