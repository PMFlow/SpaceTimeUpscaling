%% Contours - solutions of the Monod model
function plot_surf_c(x,z,c1,c2,N)

c1=c1/N; c2=c2/N;
levels=24;

figure;
contourf(x,z,c1,levels,'EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$c_1(x,y,t)$','Interpreter','latex'); 

figure;
contourf(x,z,c2,levels,'EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$c_2(x,y,t)$','Interpreter','latex'); 
