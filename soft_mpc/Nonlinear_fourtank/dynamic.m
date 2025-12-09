function x_new=dynamic(x,u,param)
%RK4 dynamics
k1 = fun_c(x,u,param);
k2 = fun_c(x+param.delta/2*k1,u,param);
k3 = fun_c(x+param.delta/2*k2,u,param);
k4 = fun_c(x+param.delta*k3,u,param);
x_new = x+param.delta/6*(k1+2*k2+2*k3+k4);
end