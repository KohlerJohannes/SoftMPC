%Generates random disturbances, with bound param.w_max and time-dependended
%scaling, as detalied in the paper
%%
param=Tank_param;
Tsim=param.t_switch(end);
%gaussian
%w=mvnrnd(zeros(param.n,1),param.Sigma,ceil(Tsim/param.delta))';
%uniform
w=(rand(param.n,ceil(Tsim/param.delta))-0.5)*2*param.w_max;
%set to 0 before t_disturbance
t=0:param.delta:Tsim;
for k=1:length(t)%Tsim/param.delta
   if t(k)<=param.disturbances_t(1) ||param.delta*k>=param.disturbances_t(3)
      w(:,k)=zeros(param.n,1);
   elseif t(k)<=param.disturbances_t(2)
       w(:,k)=w(:,k)*(t(k)-param.disturbances_t(1))/(param.disturbances_t(2)-param.disturbances_t(1));
   elseif t(k)<=param.disturbances_t(3)
       w(:,k)=w(:,k)*(param.disturbances_t(3)-t(k))/(param.disturbances_t(3)-param.disturbances_t(2));
   else 
       error('?')
   end
end
%uniform distribution
save('Results/w','w')