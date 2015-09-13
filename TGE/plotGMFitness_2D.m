close all, clear all
[FileName,PathName,FilterIndex] = uigetfile('*.txt');
A = dlmread([PathName FileName]);
NV = 2;
B = A(:,1:2*NV);
U = A(:,5:end);
solution = input('Enter with the known solution: '); % example: @(x,t) cos(x-t)

limx = [min(min(B(:,[1,3]))) max(max(B(:,[1,3])))];
limt = [min(min(B(:,[2,4]))) max(max(B(:,[2,4])))];




subplot(3,3,[1,2,3,4,5,6])
npt = 2500;
npt = ceil(sqrt(npt));
[x,t] = meshgrid(linspace(limx(1),limx(2),npt),linspace(limt(1),limt(2),npt));
mesh(x,t,solution(x,t))
colormap('gray')

hold on
for d = 1:size(U,1)
%     if (d ~= 7 && d ~= 8 && d ~= 9)
        map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
        map_xi_t = @(t) 2.*(t-B(d,3))./(B(d,4)-B(d,3)) - 1;
        [ax_x,ax_t] = meshgrid(linspace(B(d,1),B(d,2),10),linspace(B(d,3),B(d,4),10));
        mapx{1} = map_xi_x(ax_x);
        mapx{2} = map_xi_t(ax_t);
        surf(ax_x,ax_t,u_hat(U(d,:),mapx),'LineWidth',2);
%     end
end
% set(gca,'YLim',[0.1 0.8]);
hold off
%axis equal
subplot(3,3,7)
npt = 2500;
npt = ceil(sqrt(npt));
[x,t] = meshgrid(linspace(limx(1),limx(2),npt),linspace(limt(1),limt(2),npt));
mesh(x,t,solution(x,t))
colormap('gray')

hold on
for d = 2:2
%     if (d ~= 7 && d ~= 8 && d ~= 9)
        map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
        map_xi_t = @(t) 2.*(t-B(d,3))./(B(d,4)-B(d,3)) - 1;
        [ax_x,ax_t] = meshgrid(linspace(B(d,1),B(d,2),10),linspace(B(d,3),B(d,4),10));
        mapx{1} = map_xi_x(ax_x);
        mapx{2} = map_xi_t(ax_t);
        surf(ax_x,ax_t,u_hat(U(d,:),mapx),'LineWidth',2);
%     end
end
% set(gca,'YLim',[0.1 0.8]);
hold off


subplot(3,3,8)
npt = 2500;
npt = ceil(sqrt(npt));
[x,t] = meshgrid(linspace(limx(1),limx(2),npt),linspace(limt(1),limt(2),npt));
mesh(x,t,solution(x,t))
colormap('gray')

hold on
for d = 5:5
%     if (d ~= 7 && d ~= 8 && d ~= 9)
        map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
        map_xi_t = @(t) 2.*(t-B(d,3))./(B(d,4)-B(d,3)) - 1;
        [ax_x,ax_t] = meshgrid(linspace(B(d,1),B(d,2),10),linspace(B(d,3),B(d,4),10));
        mapx{1} = map_xi_x(ax_x);
        mapx{2} = map_xi_t(ax_t);
        surf(ax_x,ax_t,u_hat(U(d,:),mapx),'LineWidth',2);
%     end
end
% set(gca,'YLim',[0.1 0.8]);
hold off

subplot(3,3,9)
npt = 2500;
npt = ceil(sqrt(npt));
[x,t] = meshgrid(linspace(limx(1),limx(2),npt),linspace(limt(1),limt(2),npt));
mesh(x,t,solution(x,t))
colormap('gray')

hold on
for d = 8:8
%     if (d ~= 7 && d ~= 8 && d ~= 9)
        map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
        map_xi_t = @(t) 2.*(t-B(d,3))./(B(d,4)-B(d,3)) - 1;
        [ax_x,ax_t] = meshgrid(linspace(B(d,1),B(d,2),10),linspace(B(d,3),B(d,4),10));
        mapx{1} = map_xi_x(ax_x);
        mapx{2} = map_xi_t(ax_t);
        surf(ax_x,ax_t,u_hat(U(d,:),mapx),'LineWidth',2);
%     end
end
% set(gca,'YLim',[0.1 0.8]);
hold off






% figure
% for KKK = 1:16
%     mesh(x,t,solution(x,t))
%     colormap('gray')
% 
%     hold on
%     for d = 1:size(U,1)
%         if (d ~= KKK)
%             map_xi_x = @(x) 2.*(x-B(d,1))./(B(d,2)-B(d,1)) - 1;
%             map_xi_t = @(t) 2.*(t-B(d,3))./(B(d,4)-B(d,3)) - 1;
%             [ax_x,ax_t] = meshgrid(linspace(B(d,1),B(d,2),10),linspace(B(d,3),B(d,4),10));
%             mapx{1} = map_xi_x(ax_x);
%             mapx{2} = map_xi_t(ax_t);
%             surf(ax_x,ax_t,u_hat(U(d,:),mapx),'LineWidth',2);
%         end
%     end
%     title(KKK)
%     hold off
%     pause(2);
% end