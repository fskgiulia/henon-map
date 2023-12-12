%% Henon map
clear all
close all
x(1)=0;
y(1)=0;
a=1.25;
b=0.3;
% 10000 iterations
for i=2:10000
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
plot(x,y,'k.','MarkerSize',10)
grid on
title('Hénon map')
%% example of few iterations
clear all
close all
x(1)=3;
y(1)=0;
a=1.4;
b=0.3;
% 10000 iterations
for i=2:5
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
x
y
%% self similarity
clear all
close all
x(1)=0;
y(1)=0;
a=1.4;
b=0.3;
% 10000000 iterations
for i=2:10000000
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
figure
title('Henon map')
subplot(1,4,1)
plot(x,y,'k.','MarkerSize',4)
xlim([0.6 0.75]);
ylim([0.15 0.2]);
grid on
subplot(1,4,2)
plot(x,y,'k.','MarkerSize',4)
grid on
xlim([0.66 0.68]);
ylim([0.179 0.185]);
subplot(1,4,3)
plot(x,y,'k.','MarkerSize',4)
grid on
xlim([0.665 0.6725]);
ylim([0.1826 0.1834]);
subplot(1,4,4)
plot(x,y,'k.','MarkerSize',4)
grid on
ylim([0.183 0.1832]);
xlim([0.6705 0.6715]);
%% 10000 iterations of a point, plot the last 9000
clear all
close all
N=10000;
lim=9000;
x(1)=0;
x1(1)=0.001;
y(1)=0;
y1(1)=0.001;
for i=2:N
    x(i)=1+y(i-1)-1.4*x(i-1)^2;
    y(i)=0.3*x(i-1);
    x1(i)=1+y1(i-1)-1.4*x1(i-1)^2;
    y1(i)=0.3*x1(i-1);
end
figure
title('Last 9000 iterations of 10000')
subplot(1,3,1)
plot(x(lim:N),y(lim:N),'.r')
grid on
legend('x_0=0, y_0=0')
subplot(1,3,2)
plot(x1(lim:N),y1(lim:N),'.b')
grid on
legend('x_0=0.001, y_0=0.001')
subplot(1,3,3)
plot(x(lim:N),y(lim:N),'.r')
hold on
plot(x1(lim:N),y1(lim:N),'.b')
grid on
legend('x_0=0, y_0=0','x_0=0.001, y_0=0.001')
%% sensitivity on initial conditions
clear all
close all
N=50;
x(1)=0;
x1(1)=0.0001;
y(1)=0;
y1(1)=0.0001;
for i=2:N
    x(i)=1+y(i-1)-1.4*x(i-1)^2;
    y(i)=0.3*x(i-1);
    x1(i)=1+y1(i-1)-1.4*x1(i-1)^2;
    y1(i)=0.3*x1(i-1);
end
plot(1:N,x-x1,'.-k')
xlabel('k')
ylabel("x_k-x'_k")
grid on
%% with more iterations
clear all
close all
N=200;
x(1)=0;
x1(1)=0.0001;
y(1)=0;
y1(1)=0.0001;
for i=2:N
    x(i)=1+y(i-1)-1.4*x(i-1)^2;
    y(i)=0.3*x(i-1);
    x1(i)=1+y1(i-1)-1.4*x1(i-1)^2;
    y1(i)=0.3*x1(i-1);
end
plot(1:N,y-y1,'.-k')
grid on
%% system dynamics with varying a
clear all
close all
N=1000000;
x=zeros(1,N);
y=zeros(1,N);
a=1.4
b=0.3;
step=0.001;
amin=0.1;
figure
%subplot(1,2,1)
for i=amin:step:a
    x(1,1)=0.1;
    y(1,1)=0;
    for j=1:N
        x(1,j+1)=1-(i*x(1,j)^2)+y(1,j);
        y(1,j+1)=b*x(1,j);
    end
    plot(i*ones(1,10),x(1,N-9:N),'.k','MarkerSize',2)
    line1=sprintf('Bifurcation diagram for the Henon map with %0.f iterations',N);
    line2=sprintf('b = %.2f',b);
    title({line1,line2});
    xlabel('a');
    ylabel('x');
    hold on
    grid on
end
figure
%subplot(1,2,2)
for i=amin:step:a
    x(1,1)=0.1;
    y(1,1)=0;
    for j=1:N
        x(1,j+1)=1-(i*x(1,j)^2)+y(1,j);
        y(1,j+1)=b*x(1,j);
    end
    plot(i*ones(1,10),x(1,N-9:N),'.k','MarkerSize',2)
    line1=sprintf('Bifurcation diagram for the Henon map with %0.f iterations',N);
    line2=sprintf('b = %.2f',b);
    title({line1,line2});
    xlim([1 1.5])
    xlabel('a');
    ylabel('x');
    hold on
    grid on
end
%% system dynamics with varying b
clear all
close all
bmin=-0.6;
a=1.4;
b=0.3;
step=0.001;
N=1000000;
figure
subplot(1,2,1)
for i=bmin:step:(b+0.1)
    x(1,1)=0.1;
    y(1,1)=0;
    for j=1:N
        x(1,j+1)=1-(a*x(1,j)^2)+y(1,j);
        y(1,j+1)=i*x(1,j);
    end
    plot(i*ones(1,10),x(1,N-9:N),'.k','MarkerSize',2)
 %   line1=sprintf('Bifurcation diagram for the Henon map with %0.f iterations',N);
 %   line2=sprintf('a = %.2f',a);
 %   title({line1,line2});
    line1=sprintf('a = %.1f',a);
    title({line1})
    xlabel('b');
    ylabel('x');
    hold on
    grid on
end
subplot(1,2,2)
for i=bmin:step:(b+0.1)
    x(1,1)=0.1;
    y(1,1)=0;
    for j=1:N
        x(1,j+1)=1-(a*x(1,j)^2)+y(1,j);
        y(1,j+1)=i*x(1,j);
    end
    plot(i*ones(1,10),x(1,N-9:N),'.k','MarkerSize',2)
    xlim([-0.1 0.4])
  %  line1=sprintf('Bifurcation diagram for the Henon map with %0.f iterations',N);
  %  line2=sprintf('a = %.2f',a);
  %  title({line1,line2});
    line1=sprintf('a = %.1f',a);
    title({line1})
    xlabel('b');
    ylabel('x');
    hold on
    grid on
end
%% Henon map with different a,b
clear all 
close all
x(1)=0;
y(1)=0;
a=[-0.2, 0.2, -0.5, -0.5, 0.2, -0.22];
b=[-0.9991, 0.99912, -0.9999, -0.9995, 0.9991, -0.9999];
% 10000 iterations
k=1;
for j=1:length(a)
    for i=2:10000
        x(i)=1-a(j)*(x(i-1)^2)+y(i-1);
        y(i)=b(j)*x(i-1);
    end
    subplot(3,2,k)
    plot(x,y,'k.','MarkerSize',2)
    grid on
    line1=sprintf('Henon map');
    line2=sprintf('a = %.4f, b = %.4f', a(j),b(j));
    title({line1,line2})
    k=k+1;
end
%% Lozi
clear all
close all
x(1)=0;
y(1)=0;
a=1.7;
b=0.5;
% 10000 iterations
for i=2:100000
    x(i)=1-a*(abs(x(i-1)))+y(i-1);
    y(i)=b*x(i-1);
end
figure
plot(x,y,'k.','MarkerSize',2)
grid on
line1=sprintf('Lozi Strange Attractor');
line2=sprintf('a = %.1f, b = %.1f',a,b);
title({line1,line2})
%% enlargements
clear all
close all
x(1)=0;
y(1)=0;
a=1.4;
b=0.3;
for i=2:100000
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
figure
subplot(1,2,1)
plot(x,y,'k.','MarkerSize',2)
grid on
xlim([0 1])
ylim([0 0.3])
subplot(1,2,2)
plot(x,y,'k.','MarkerSize',2)
grid on
xlim([0.7 0.8])
ylim([0.15 0.18])
%% fixed points 
clear all
close all
a=1.4;
b=0.3;
x1=(-(1-b)+sqrt((1-b)^2+4*a))/(2*a)
y1=b*x1
x2=(-(1-b)-sqrt((1-b)^2+4*a))/(2*a)
y2=b*x2
x(1)=0;
y(1)=0;
N=10000;
for i=2:N
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
figure
subplot(1,2,1)
plot(x(1:N-1),x(2:N),'k.','MarkerSize',2)
xlabel('x_n')
ylabel('x_{n+1}')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
hold on
%axis equal
plot(x1,y1,'*r',LineWidth=2)
plot(x2,y2,'*r',LineWidth=2)
grid on
line1=sprintf('Fixed points of the Hénon map');
line2=sprintf('a = %.1f, b = %.1f', a,b);
title({line1,line2})
subplot(1,2,2)
plot(x(1:N-1),x(2:N),'g.','MarkerSize',2)
hold on
plot(x(2:N),x(1:N-1),'b.','MarkerSize',2)
xlabel('x_n')
ylabel('x_{n+1}')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
%axis equal
%plot(linspace(-2,1,2),linspace(-2,1,2),'k--',MarkerSize=1)
plot(x1,y1,'*r',LineWidth=2)
plot(x2,y2,'*r',LineWidth=2)
legend('Backward iterate','Forward iterate')
grid on
line1=sprintf('Fixed points of the Hénon map');
line2=sprintf('a = %.1f, b = %.1f', a,b);
title({line1,line2})
%% animation of the dynamics of Henon map
clear all
close all
a=1.4;
b=0.3;
x(1)=0;
y(1)=0;
figure
plot(x(1),y(1),'.',MarkerSize=8)
line1=sprintf('Hénon map')
line2=sprintf('a = %.1f, b = %.1f', a,b)
title({line1,line2})
hold on
for i=2:1000
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
    plot(x(i),y(i),'.','MarkerSize',8)
    grid on
    xlim([-1.5 1.6])
    ylim([-0.5 0.5])
    %pause(2*exp(-2*i))
end
%% fixed points for different a
clear all
close all
a=1.2;
b=0.3;
x1=(-(1-b)+sqrt((1-b)^2+4*a))/(2*a)
y1=b*x1
x2=(-(1-b)-sqrt((1-b)^2+4*a))/(2*a)
y2=b*x2
x(1)=0;
y(1)=0;
N=10000;
for i=2:N
    x(i)=1-a*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
figure
subplot(1,2,1)
plot(x(1:N-1),x(2:N),'k.','MarkerSize',2)
xlabel('x_n')
ylabel('x_{n+1}')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
hold on
%axis equal
plot(x1,y1,'*r',LineWidth=2)
plot(x2,y2,'*r',LineWidth=2)
grid on
line1=sprintf('Fixed points of the Hénon map');
line2=sprintf('a = %.1f, b = %.1f', a,b);
title({line1,line2})
subplot(1,2,2)
plot(x(1:N-1),x(2:N),'g.','MarkerSize',2)
hold on
plot(x(2:N),x(1:N-1),'b.','MarkerSize',2)
xlabel('x_n')
ylabel('x_{n+1}')
%xlim([-1.5 1.5])
%ylim([-1.5 1.5])
%axis equal
%plot(linspace(-2,1,2),linspace(-2,1,2),'k--',MarkerSize=1)
plot(x1,y1,'*r',LineWidth=2)
plot(x2,y2,'*r',LineWidth=2)
legend('Backward iterate','Forward iterate')
grid on
line1=sprintf('Fixed points of the Hénon map');
line2=sprintf('a = %.1f, b = %.1f', a,b);
title({line1,line2})
%% Henon map with different a,b part2: DIVERGENCE OF THE MAP
clear all 
close all
x(1)=0;
y(1)=0;
a=[-0.6,-0.5,-0.4,-0.3250,-0.3,-0.2, -0.15, 0, 0.05, 0.2];
c=[-0.7,-0.5,-0.999,-0.999,-0.999,0,0,0.2,0.3,0.6];
% 10000 iterations
for j=1:length(a)
    figure
    b=linspace(-0.9991,c(j),100);
    for k=1:length(b)
        for i=2:10000
            x(i)=1-a(j)*(x(i-1)^2)+y(i-1);
            y(i)=b(k)*x(i-1);
        end
        plot(x,y,'k.','MarkerSize',2)
        xlabel('x')
        ylabel('y')
        grid on
        line1=sprintf('Henon map');
        line2=sprintf('a = %.4f, b = [-0.9991, %.4f]', a(j),b(k));
        title({line1,line2})
        hold on
        pause(0.01)
    end
   % pause(2)
  %  close all
end
%% Henon map with different a,b part3: DIVERGENCE OF THE MAP
clear all 
close all
x(1)=0;
y(1)=0;
a=[-0.67,-0.66,-0.63,-0.57,-0.1999, -0.001, 0.0001,0.13, 0.2,0.5,0.9,1];
c=[-0.99,-0.9,-0.9,-0.9,0,0,0.2,0.4,0.6,0.8,0.5,0.5,0.45];
% 10000 iterations
for j=1:length(a)
    figure
    b=linspace(-0.9999,c(j),100);
    for k=1:length(b)
        for i=2:10000
            x(i)=1-a(j)*(x(i-1)^2)+y(i-1);
            y(i)=b(k)*x(i-1);
        end
        plot(x,y,'k.','MarkerSize',2)
        xlabel('x')
        ylabel('y')
        grid on
        line1=sprintf('Henon map');
        line2=sprintf('a = %.4f, b = [-0.9999, %.4f]', a(j),b(k));
        title({line1,line2})
        hold on
        pause(0.01)
    end
  %  pause(2)
  %  close all
end
%% eigenvalues - fixed points
clear all
close all
H1=[2.8*0.6314, 1; 0.3, 0]
eig(H1)
det(H1)
H2=[2.8*-1.1314, 1; 0.3, 0]
eig(H2)