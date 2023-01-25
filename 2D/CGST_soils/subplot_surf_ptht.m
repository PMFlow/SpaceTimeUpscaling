%% Contours of flow solutions

function subplot_surf_ptht(x,z,p,tht,Vx,Vy)

levels=24;

%% one plot
figure;

subplot(1,3,1)
contourf(x,z,p,levels,'EdgeColor','none'); colorbar(ticks=[-1 -0.8 -0.6 -0.4 -0.2 0]); colormap(parula); view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$\psi(x,z,t)$','Interpreter','latex'); 

subplot(1,3,2)
contourf(x,z,tht,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$\theta(x,z,t)$','Interpreter','latex'); 


subplot(1,3,3)
contourf(x,z,sqrt(Vx.^2+Vy.^2),levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$q(x,z,t)$','Interpreter','latex'); 


%% 4 plots
% figure;
% contourf(x,z,p,levels,'EdgeColor','none'); colorbar(ticks=[-1 -0.8 -0.6 -0.4 -0.2 0]); colormap(parula); view(0,90)
% set(gca, 'XTick', []); set(gca, 'YTick', []);
% title('$\psi(x,z,t)$','Interpreter','latex'); 
% 
% figure;
% contourf(x,z,tht,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
% set(gca, 'XTick', []); set(gca, 'YTick', []);
% title('$\theta(x,z,t)$','Interpreter','latex'); 
% 
% figure;
% contourf(x,z,Vx,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
% set(gca, 'XTick', []); set(gca, 'YTick', []);
% title('$q_x(x,z,t)$','Interpreter','latex'); 

% figure;
% contourf(x,z,Vy,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
% set(gca, 'XTick', []); set(gca, 'YTick', []);
% title('$q_z(x,z,t)$','Interpreter','latex'); 
