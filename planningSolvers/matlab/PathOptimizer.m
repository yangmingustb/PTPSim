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
s1 = s0+reso_s*1;
s2 = s0+reso_s*2;

%{
xr(s1)=s1;
yr(s1)=0;
thetar(s1)=0;

xr(s2)=s2;
yr(s2)=0;
thetar(s2)=0;
%}

rho1 = (tan(theta0)*xr(s1)-tan(theta0)*x0+y0-yr(s1))/(sin(thetar(s1)+pi/2)-tan(theta0)*cos(thetar(s1)+pi/2));

x1 = xr(s1) + rho1*cos(thetar(s1)+pi/2);
y1 = yr(s1) + rho1*sin(thetar(s1)+pi/2);
m = kappa0*((x1-x0)^2+(y1-y0)^2)^(3/2);
n = (x1-x0)*(yr(s2)-2*y1+y0)+(y1-y0)*(xr(s2)-2*x1+x0);
q = (x1-x0)*sin(thetar(s2)+pi/2)-(y1-y0)*cos(thetar(s2)+pi/2);
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


H1 = eye(n_s+1);
i=1;
while(i<n_s+1)
    H1(i,i+1)=-1;
    H1(i+1,i)=-1; 
    i = i+1;
end
tmp = zeros(n_s+1,1);
H1 = [H1,tmp];
tmp = zeros(1, n_s+2);
H1 = [H1;tmp];

H2 = eye(n_s+2);
i=1;
while(i<=n_s+2)
    if i ==1
        
        H2(i,i+1)=-2;
        H2(i,i+2)=1;
        H2(i+1,i)=-2; 
        H2(i+2,i)=1; 
    end
    
    if i == 2
        H2(i,i)=5;
        H2(i,i+1)=-4;
        H2(i,i+2)=1;
        H2(i+1,i)=-4; 
        H2(i+2,i)=1; 
    end
    
    if i >=3 || i <=n_s
        H2(i,i)=6;
        H2(i,i+1)=-4;
        H2(i,i+2)=1;
        H2(i+1,i)=-4; 
        H2(i+2,i)=1; 
    end
    if i == n_s+1
        H2(i,i)=5;
        H2(i,i+1)=-2;
        H2(i+1,i)=-2;
    end     
    i = i+1;
end

H3 = eye(n_s);
tmp = zeros(n_s,2);
H3 = [H3,tmp];
tmp = zeros(2, n_s+2);
H3 = [H3;tmp];

H = 2*w_smt/reso_s*H1+2*ww_smt/(reso_s^3)*H2+2*reso_s*w_ref*H3;
q = -2*reso_s * w_ref*rho_r;

i=1;
rho = [];
while(i<=n_s+2)
    
    rho = [rho; rho(i)];
    i=i+1;
end

f = 1/2*rho'.*H.*rho+q.*rho;

    
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
