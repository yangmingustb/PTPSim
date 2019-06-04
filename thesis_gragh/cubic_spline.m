start=[0,0];
m=[4,8,12,16];
n=[-2,-1,0,1,2];

%第一组
i=1;
while i<=length(n)
    
    x=[start(1),m(1)];
    y=[start(2),n(i)];
    plot_graph(x,y);
    i=i+1;
end

j=1;
while j<=length(n)
    q=1;
    while q<=length(n)
    
        x=[m(1),m(2)];
        y=[n(j),n(q)];
        plot_graph(x,y);
        q=q+1;
    end
    j=j+1;
end

%第三组
j=1;
while j<=length(n)
    q=1;
    while q<=length(n)
    
        x=[m(2),m(3)];
        y=[n(j),n(q)];
        plot_graph(x,y);
        q=q+1;
    end
    j=j+1;
end


%第四组
j=1;
while j<=length(n)
    q=1;
    while q<=length(n)
    
        x=[m(3),m(4)];
        y=[n(j),n(q)];
        plot_graph(x,y);
        q=q+1;
    end
    j=j+1;
end
xlabel('longitudinal station');
ylabel('lateral offset');

%画出单条曲线
function f=plot_graph(x,y)
    %x=[0 2];
    %y=[0 1];
    cs=spline(x,[0 y 0]);
    xx=linspace(x(1),x(2),40);
    plot(x,y,'or',xx,ppval(cs,xx),'-b','LineWidth',0.6);
    axis('equal'); 
    hold on;
end