%%   2019-4-24-15:03
%{
min f(x)
s.t. A x<=B
Aeq x=Beq
C(x)<=0
Ceq(x)=0
LB<=X<=UB


%}
clc;

show_animation = 1;

%%  parameters
global rho_max;
rho_max = 3.0;    %   ��ʻͨ�������ұ߽�
global rho_min;
rho_min = -3.0;

global s_max;
s_max = 100.0;
global reso_s;
reso_s = 1.0;
global n_s;
n_s = s_max/reso_s;
disp(n_s);
global r_circle;
r_circle = 1.0;   %   �������Բ
global d_circle;
d_circle = 2.0;   %   Բ�ļ��
global obs_inflation;    %   �ϰ�������
obs_inflation = 1.0;
global kappa_max;
kappa_max = 0.187;
global w_d;
w_d = 1.0;
global w_dd;
w_dd = 10.0;
global w_ddd;
w_ddd = 50.0;
global w_ref;
w_ref = 1.5;
global rho_r;
rho_r = zeros(1, n_s+2);

%%  ���״̬
x0=0;
y0=0;
s0=0;
rho0=0;
theta0=0;
kappa0=0;
x_init = [x0,y0,s0,rho0,theta0,kappa0];

%%  ���s������

global s;
i=1;
s = reso_s:reso_s:(s_max+3*reso_s);
%   disp(s);

%%  �ϰ����ʾ
global static_obs;
static_obs = [60, 30; 1,-1];    %   �ϰ����frenet����,
[r,c] = size(static_obs);
global n_static;
n_static = c;   %   ��̬�ϰ������

%%  ��ʽԼ��
rho1 = reso_s*(tan(theta0-thetar(s0)))+rho0;

x1 = xr(s(1)) + rho1*cos(thetar(s(1))+pi/2);
y1 = yr(s(1)) + rho1*sin(thetar(s(1))+pi/2);
m = kappa0*((x1-x0)^2+(y1-y0)^2)^(3/2);
n = (x1-x0)*(yr(s(2))-2*y1+y0)+(y1-y0)*(xr(s(2))-2*x1+x0);
q = (x1-x0)*sin(thetar(s(2))+pi/2)-(y1-y0)*cos(thetar(s(2))+pi/2);
rho2 = (m-n)/q;
%   disp(rho1);
%   disp(rho2);

%%  SQP��������
rho_init_guess = zeros(n_s+3, 1);
A = [];
b = [];   % A*x <= b
Aeq = zeros(2,n_s+3);
Aeq(1,1)=1;
Aeq(2,2)=1;
Beq = [rho1; rho2];  % Aeq*x = Beq
lb = rho_min * ones(1, n_s+3);
ub = rho_max * ones(1, n_s+3);
%   disp(lb);
%   disp(ub);

%%  �㷨����
%   interior-point,trust-region-reflective,
%   sqp,sqp-legacy,active-et
%   �ڵ㷨���
tic
options = optimoptions('fmincon','Algorithm', 'interior-point');
%options = optimoptions('fmincon','Algorithm', 'sqp');
%options = optimoptions('fmincon','Algorithm', 'sqp-legacy');
%options = optimoptions('fmincon','Algorithm', 'active-set');    %����

[rho_sqp,f_sqp, exitflag, output]=fmincon(@obj,rho_init_guess,A,b,Aeq,Beq,lb,ub,@nonlcon,options);
sqp_time = toc;
disp("sqp_time:");
disp(sqp_time);
%   disp("SQP Result: ");
%   disp(rho_sqp);
size_rho = size(rho_sqp);
disp(size_rho);

%%  plot gragh
if(show_animation)
    figure(1);
    clf;
    hold on;
    subplot(211);
    hold on;
    s = [s0,s];
    rho_sqp = [rho0, rho_sqp'];
    plot(s,rho_sqp);
    xlabel('station');
    ylabel('offset');
    axis([0, 110, -6, 6]);
    hold on;

    scatter(static_obs(1,1),static_obs(2,1),30,'r');
    scatter(static_obs(1,2),static_obs(2,2),30,'r');
    
    axis("equal");
    tmp_kappa = trajectory_kappa(rho_sqp);
    max_kappa = max(tmp_kappa);
    min_kappa = min(tmp_kappa);
    disp("max_kappa:");
    disp(max_kappa);
    disp("min_kappa:");
    disp(min_kappa);
    subplot(212);
    hold on;
    tmp_kappa = [kappa0, tmp_kappa];
    %plot(s(1:(n_s+1)), tmp_kappa);
    hold on;

    %% ����ͼƽ��������Ӧ�ü���
    %
    si = s0:0.01:s_max;
    yi=interp1(s(1:(n_s+1)),tmp_kappa,si,'spline');
    plot(si,yi);
    axis([0, 110, -1, 1]);
    xlabel('station');
    ylabel('curvature');
        
    %}
    hold on;
    k_max = kappa_max*ones(1, n_s+4);
    plot(s,k_max,'-.r','LineWidth',1);
    plot(s,-k_max,'-.r','LineWidth',1);
end

%%  Objective function
function f=obj(rho)
global reso_s;
global w_d;
global w_dd;
global w_ddd;
global w_ref;
global rho_r;
global n_s;

f_d = 0;
f_dd = 0;
f_ddd = 0;
f_ref = 0;
i=1;
while(i<=n_s)
    f_d = f_d + (rho(i+1)-rho(i))^2/reso_s^2;
    f_dd = f_dd +(rho(i+2)-2*rho(i+1)+rho(i))^2/reso_s^4;
    f_ddd = f_ddd + (rho(i+3)-3*rho(i+2)+3*rho(i+1)-rho(i))^2/reso_s^6;
    f_ref = f_ref + (rho(i)-rho_r(i))^2;
    i= i+1;
end
%{
disp('=========================================');
disp(f_d);
disp(f_dd);
disp(f_ddd);
disp(f_ref);
%}

f_d = reso_s * w_d * f_d;
f_dd = reso_s * w_dd * f_dd;
f_ddd = reso_s * w_ddd * f_ddd;
f_ref = reso_s * w_ref * f_ref;

f = f_d + f_dd + f_ddd + f_ref;

end

%%  Nonlinear Constraint Function 
function [c, ceq] = nonlcon(rho)
global kappa_max;
global s;
global n_s;
global r_circle;
global d_circle;
global static_obs;
global obs_inflation;
c_kappa1 = zeros(n_s, 1);
c_kappa2 = zeros(n_s, 1);
c1_obs = zeros(n_s+3, 1);
c2_obs = zeros(n_s+3, 1);
c3_obs = zeros(n_s+3, 1);

%%  ����Լ��
i=1;
while(i<=n_s)
    
    k1 = (x(s(i+1),rho(i+1))-x(s(i),rho(i)))*(y(s(i+2),rho(i+2))-2*y(s(i+1),rho(i+1))+y(s(i),rho(i)));
    k2 = (y(s(i+1),rho(i+1))-y(s(i),rho(i)))*(x(s(i+2),rho(i+2))-2*x(s(i+1),rho(i+1))+x(s(i),rho(i)));
    k3 = ((x(s(i+1),rho(i+1))-x(s(i),rho(i)))^2+(y(s(i+1),rho(i+1))-y(s(i),rho(i)))^2)^(3/2);
    c_kappa1(i) = k1-k2-k3*kappa_max;
    c_kappa2(i) = -k1+k2-k3*kappa_max; 
    i=i+1;
    
end

%%  �ϰ���Լ��
i=1;
while(i<=(n_s+3))
    c1_obs(i) = (r_circle+obs_inflation)^2-(s(i)-static_obs(1,1))^2-(rho(i)-static_obs(2,1))^2;
    c2_obs(i) = (r_circle+obs_inflation)^2-(s(i)+d_circle-static_obs(1,1))^2-(rho(i)-static_obs(2,1))^2;
    c3_obs(i) = (r_circle+obs_inflation)^2-(s(i)+2*d_circle-static_obs(1,1))^2-(rho(i)-static_obs(2,1))^2;
    i=i+1;
end

c1_obs2 = zeros(n_s+3, 1);
c2_obs2 = zeros(n_s+3, 1);
c3_obs2 = zeros(n_s+3, 1);
i=1;
while(i<=(n_s+3))
    c1_obs2(i) = (r_circle+obs_inflation)^2-(s(i)-static_obs(1,2))^2-(rho(i)-static_obs(2,2))^2;
    c2_obs2(i) = (r_circle+obs_inflation)^2-(s(i)+d_circle-static_obs(1,2))^2-(rho(i)-static_obs(2,2))^2;
    c3_obs2(i) = (r_circle+obs_inflation)^2-(s(i)+2*d_circle-static_obs(1,2))^2-(rho(i)-static_obs(2,2))^2;
    i=i+1;
end

c = [c_kappa1;c_kappa2;c1_obs;c2_obs;c3_obs;c1_obs2;c2_obs2;c3_obs2];
ceq = [];
end

%%  FrenetToCartesian(s, rho)
function [tmp_x] = x(s,rho)
tmp_x = xr(s) + rho*cos(thetar(s)+pi/2);
end
function [tmp_y] = y(s,rho)
tmp_y = yr(s) + rho*sin(thetar(s)+pi/2);
end

%%  ������ֱ��
function [tmp_xr] = xr(s)
tmp_xr = s;
end
function [tmp_yr] = yr(s)
tmp_yr = 0*s;
end
function [tmp_thetar] = thetar(s)
tmp_thetar = 0*s;
end

%%  kappa���㺯��
function [tmp_kappa]=trajectory_kappa(rho)
global n_s;
global s;
tmp_kappa = zeros(1, n_s);
i=1;
while(i<=n_s)
    k1 = (x(s(i+1),rho(i+1))-x(s(i),rho(i)))*(y(s(i+2),rho(i+2))-2*y(s(i+1),rho(i+1))+y(s(i),rho(i)));
    k2 = (y(s(i+1),rho(i+1))-y(s(i),rho(i)))*(x(s(i+2),rho(i+2))-2*x(s(i+1),rho(i+1))+x(s(i),rho(i)));
    k3 = ((x(s(i+1),rho(i+1))-x(s(i),rho(i)))^2+(y(s(i+1),rho(i+1))-y(s(i),rho(i)))^2)^(3/2);
    if(k3 == 0.0)
        tmp_kappa(i) = 0;
    
    else
        tmp_kappa(i) = (k1-k2)/k3;
    end
    i=i+1;

end
end

