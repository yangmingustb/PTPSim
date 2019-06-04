x=1:1:13;
y=[0.5,0.5,0.5,1,1,0,1,1,1,0.2,0.2,0.2,0.2];
disp(y(1));
stem(x,y,"fill",'b','LineWidth',1);
axis([0, 14, 0, 1.4]);
set(gca,'XTick',0:1:13);
set(gca,'YTick',[0 0.2 0.5 1]);
xlabel('path index');
ylabel('static collision check');
%set(gca,'YTickLabel',{'0.2','0.5','1'});
%set(gca,'Fontname','Times newman','FontSize',10);
%set(gcf,'Units','centimeters','Position',[10 5 7 7]);
%set(get(gca,'XLabel'),'FontSize',10);
%set(get(gca,'YLabel'),'FontSize',10);
%set(gca,'Position',[.15 .15 .8 .75]);


