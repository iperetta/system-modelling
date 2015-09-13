close all, clear all
[FileName,PathName,FilterIndex] = uigetfile('*.txt');
A = dlmread([PathName FileName]);
NV = 1;
B = A(:,1:2*NV);
U = A(:,3:end);

solution = input('Enter with the known solution: '); % example: @(x) 1./x.*erf(sqrt(2)/2*x)
ma = min(B(:)); mb = max(B(:));
ax_x = linspace(ma,mb,100);
subplot(3,3,[1,2,3,4,5,6])
plot(ax_x,solution(ax_x),'-','LineWidth',4,'Color',[0.75,0.75,0.75]);
liminf = min(solution(ax_x));
limsup = max(solution(ax_x));
hold on
for d = 1:size(U,1)
    map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
    ax_x = linspace(B(d,1),B(d,2),100);
    mapx{1} = map_xi_x(ax_x);
    plot(ax_x,u_hat(U(d,:),mapx),'--k','LineWidth',2);
    set(gca,'XLim',[ma mb],'YLim',[-0.5 1]);
    scatter(B(d,1),solution(B(d,1)),'>k')
    scatter(B(d,2),solution(B(d,2)),'<k')
end
legend('Expected solution','Piecewise polynomial approximations','Piecewise domain boundary [a_i,b_i]')
xlabel('(a)','FontSize',14)
set(gca,'FontSize',14)
grid on
hold off

subplot(3,3,7)
hold on
ax_x = linspace(ma,mb,100);
plot(ax_x,solution(ax_x),'-','LineWidth',4,'Color',[0.75,0.75,0.75]);
d = 2;
map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
ax_x = linspace(B(d,1),B(d,2),100);
mapx{1} = map_xi_x(ax_x);
plot(ax_x,u_hat(U(d,:),mapx),'--k','LineWidth',2);
set(gca,'XLim',[ma mb],'YLim',[-0.5 1]);
scatter(B(d,1),solution(B(d,1)),'>k')
scatter(B(d,2),solution(B(d,2)),'<k')
xlabel('(b)','FontSize',14)
set(gca,'FontSize',14)
grid on
hold off

subplot(3,3,8)
hold on
ax_x = linspace(ma,mb,100);
plot(ax_x,solution(ax_x),'-','LineWidth',4,'Color',[0.75,0.75,0.75]);
d = 7;
map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
ax_x = linspace(B(d,1),B(d,2),100);
mapx{1} = map_xi_x(ax_x);
plot(ax_x,u_hat(U(d,:),mapx),'--k','LineWidth',2);
set(gca,'XLim',[ma mb],'YLim',[-0.5 1]);
scatter(B(d,1),solution(B(d,1)),'>k')
scatter(B(d,2),solution(B(d,2)),'<k')
xlabel('(c)','FontSize',14)
set(gca,'FontSize',14)
grid on
hold off

subplot(3,3,9)
hold on
ax_x = linspace(ma,mb,100);
plot(ax_x,solution(ax_x),'-','LineWidth',4,'Color',[0.75,0.75,0.75]);
d = 11;
map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
ax_x = linspace(B(d,1),B(d,2),100);
mapx{1} = map_xi_x(ax_x);
plot(ax_x,u_hat(U(d,:),mapx),'--k','LineWidth',2);
set(gca,'XLim',[ma mb],'YLim',[-0.5 1]);
scatter(B(d,1),solution(B(d,1)),'>k')
scatter(B(d,2),solution(B(d,2)),'<k')
xlabel('(d)','FontSize',14)
set(gca,'FontSize',14)
grid on
hold off