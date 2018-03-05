% This is to plot squared-Dirichlet density defined on sphere

clear;
addpath('~/Documents/MATLAB/tight_subplot/');

% plot for different parameters
alpha=[1,0.5,0.1];
L=length(alpha);

% get the coordinates of a unit sphere
N_pts=60;
[x,y,z]=sphere(N_pts);

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,L,[0,.07],[.11,.1],[.05,.05]);

for i=1:L
    v=(2*alpha(i)-1)*sum(log(abs([x(:),y(:),z(:)]+1e-16)),2);
    v=exp(reshape(v,N_pts+1,N_pts+1));
    subplot(ha(i));
    surf(x,y,z,v,'parent',ha(i)), shading interp
    set(gca,'fontsize',14);
    xlabel('l_1','fontsize',14);ylabel('l_2','fontsize',14);zlabel('l_3','fontsize',14,'rotation',0);
    title(['\alpha = ',num2str(alpha(i))],'fontsize',16);
end

fig.PaperPositionMode = 'auto';
print(fig,'./result/sqdirichlet','-dpng','-r0');