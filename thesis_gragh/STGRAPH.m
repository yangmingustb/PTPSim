axis([0, 10, 0, 120]);
set(gca,'xtick',0:0.5:10);
set(gca,'ytick',0:2:120);
grid on;
%set(gca,'GridLineStyle','-','GridColor','k','GridAlpha',0.5);
hold on;
p0=[0,0];
  
p1=[0.5,0];
p2=[0.5,2];
p3=[0.5,4];
p4=[0.5,6];
p5=[0.5,8];
p6=[0.5,10];


plot([p0(1),p1(1)],[p0(2),p1(2)],'r','LineWidth',1);
plot([p0(1),p2(1)],[p0(2),p2(2)],'r','LineWidth',1);
plot([p0(1),p3(1)],[p0(2),p3(2)],'r','LineWidth',1);
plot([p0(1),p4(1)],[p0(2),p4(2)],'r','LineWidth',1);
plot([p0(1),p5(1)],[p0(2),p5(2)],'r','LineWidth',1);
plot([p0(1),p6(1)],[p0(2),p6(2)],'r','LineWidth',1);


p7=[1,10];
p8=[1,12];
p9=[1,14];
p10=[1,16];
p11=[1,18];
p12=[1,20];

plot([p6(1),p7(1)],[p6(2),p7(2)],'r','LineWidth',1);
plot([p6(1),p8(1)],[p6(2),p8(2)],'r','LineWidth',1);
plot([p6(1),p9(1)],[p6(2),p9(2)],'r','LineWidth',1);
plot([p6(1),p10(1)],[p6(2),p10(2)],'r','LineWidth',1);
plot([p6(1),p11(1)],[p6(2),p11(2)],'r','LineWidth',1);
plot([p6(1),p12(1)],[p6(2),p12(2)],'r','LineWidth',1);

scatter(p1(1),p1(2),5,'bo');
scatter(p2(1),p2(2),5,'bo');
scatter(p3(1),p3(2),5,'bo');
scatter(p4(1),p4(2),5,'bo');
scatter(p5(1),p5(2),5,'bo');
scatter(p6(1),p6(2),5,'bo');
scatter(p7(1),p7(2),5,'bo');
scatter(p8(1),p8(2),5,'bo');
scatter(p9(1),p9(2),5,'bo');
scatter(p10(1),p10(2),5,'bo');
scatter(p11(1),p11(2),5,'bo');
scatter(p12(1),p12(2),5,'bo');


