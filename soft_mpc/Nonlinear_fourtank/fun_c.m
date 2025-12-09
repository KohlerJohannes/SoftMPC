
function f=fun_c(x,u,param)
%cont. time dynamics
%define state-inputs
%make sure non-negative
h_1=x(1); h_2=x(2); h_3=x(3); h_4=x(4);
q_a=u(1);q_b=u(2);
f=[-param.c_1*sqrt(h_1)+param.c_3*sqrt(h_3)+param.c_u1*q_a;...
   -param.c_2*sqrt(h_2)+param.c_4*sqrt(h_4)+param.c_u2*q_b;...
   -param.c_3*sqrt(h_3)+param.c_u3*q_b;...
   -param.c_4*sqrt(h_4)+param.c_u4*q_a];

end 