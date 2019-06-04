clc;

%%  parameters
%先假设一个直道，东西方向
rho_max = 3;    %   行驶通道的左右边界
rho_min = -3;

s_max = 100;
reso_s = 1;
n_s = int16(s_max/reso_s);
%   disp(n_s);
kappa_max = 0.187;
w_smt = 1;
w_ref = 1;
rho_r = zeros(1, n_s+2);

x0=0;
y0=0;
s0=0;
rho0=0;
theta0=0;
kappa0=0;
x_init = [x0,y0,s0,rho0,theta0,kappa0];

%%
i=1;
while(i<=n_s)
    s(i)=s0+i*reso_s;
    i=i+1;
end

%{
xr(s1)=s1;
yr(s1)=0;
thetar(s1)=0;

xr(s2)=s2;
yr(s2)=0;
thetar(s2)=0;
%}

rho1 = (tan(theta0)*xr(s(1))-tan(theta0)*x0+y0-yr(s(1)))/(sin(thetar(s(1))+pi/2)-tan(theta0)*cos(thetar(s(1))+pi/2));

x1 = xr(s(1)) + rho1*cos(thetar(s(1))+pi/2);
y1 = yr(s(1)) + rho1*sin(thetar(s(1))+pi/2);
m = kappa0*((x1-x0)^2+(y1-y0)^2)^(3/2);
n = (x1-x0)*(yr(s(2))-2*y1+y0)+(y1-y0)*(xr(s(2))-2*x1+x0);
q = (x1-x0)*sin(thetar(s(2))+pi/2)-(y1-y0)*cos(thetar(s(2))+pi/2);
rho2 = (m-n)/q;
%   disp(rho1);
%   disp(rho2);

%%
x0 = [];
A = [];
b = [];   % A*x <= b
Aeq = zeros(2,n_s+2);
Aeq(1,1)=1;
Aeq(2,2)=1;

Beq = [rho1; rho2];  % Aeq*x = Beq
lb = rho_min * ones(1, n_s+2);
ub = rho_max * ones(1, n_s+2);
%   disp(lb);
%   disp(ub);


%%
tic
options = optimoptions('fmincon','Algorithm', 'sqp');
[x_sqp,f_sqp, exitflag, output]=fmincon(@obj,x0,A,b,Aeq,Beq,lb,ub,@nonlcon,options);
sqp_time = toc;
disp("sqp_time:");
disp(sqp_time);
disp("SQP Result: ");
disp(x_sqp);

%%  Objective function
function f=obj(rho)

f_smt = 0;
f_ref = 0;
i=1;
while(i<=n_s)
    f_smt = f_smt + (rho(i+1)-rho(i))^2/reso_s^2+(rho(i+2)-2*rho(i+1)+rho(i))^2/reso_s^4;
    f_ref = f_ref + (rho(i)-rho_r(i))^2;
    i= i+1;
end
f_smt = reso_s*w_smt*f_smt;
f_ref = reso_s*w_ref*f_ref;

f = f_smt+f_ref;
    
end

%%  Nonlinear Constraint Function 
function [c, ceq] = nonlcon(rho)

c = -x1*x2-10;

ceq = [];
end

%%  FrenetToCartesian(s, rho)
function [x] = x(s,rho)
x = xr(s) + rho*cos(thetar(s)+pi/2);
end
function [y] = y(s,rho)
y = yr(s) + rho*sin(thetar(s)+pi/2);
end


%%  这里是直道
function [xr] = xr(s)
xr = s;
end
function [yr] = yr(s)
yr = 0*s;
end
function [thetar] = thetar(s)
thetar = 0*s;
end
%%  
