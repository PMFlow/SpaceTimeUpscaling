function plot_CG(strvect,xx,x,c1,mc1,mc1p,mc1v,mxx1,mvv1,mxv1,D,N)

NameArray = {'Marker'}; ValueArray = {'o','+','*'}'; 
figure
P=plot(x,mc1p/N); hold all; plot(xx,c1/N,'g-',LineWidth = 1.5);
set(P,NameArray,ValueArray,LineWidth = 1.5); 
xlabel('$x$','Interpreter','latex');
ylabel('$2a<1>/\,(N_a\mathcal{N})\,(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mc1/N); hold; plot(x,mc1v/N,'g-',LineWidth = 1.5);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$<1>/\,\mathcal{N}\,(x,t)$','Interpreter','latex'); 
set(gca, 'XTick', [0:0.2:1]);
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mxx1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{r}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mvv1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{\xi}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mxv1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{x\xi}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,(-mxv1./mc1/D-1)*100);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$(cg\_D-D)/D$ [\%]','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');
