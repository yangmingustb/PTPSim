%%  曲率与车轮转角的关系图
L = 2.7;
k = 0:0.001:0.187;
delta = atand(k*L);

figure(1);
clf;
hold on;
subplot(211);
hold on;
plot(k,delta);
xlabel('k');
ylabel('delta');
hold on;
%axis("equal");

%% speed and kappa,侧滑条件
subplot(212);
hold on;
u=0.5;
g=9.8;

v = sqrt(u*g./k);
plot(k,v);
xlabel('k');
ylabel('speed');
