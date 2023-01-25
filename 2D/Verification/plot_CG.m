function plot_CG(strvect,xx,x,c1,mc1,mc1p,mc1v,mxx1,myy1,mVx1,mVy1,mxVx1,mxVy1,myVx1,myVy1,N)

NameArray = {'Marker'}; ValueArray = {'o','+','*'}'; 
figure
P=plot(x,mc1p/N); hold all; plot(xx,c1/N,'g-',LineWidth = 1.5);
set(P,NameArray,ValueArray,LineWidth = 1.5); 
xlabel('$x$','Interpreter','latex');
ylabel('$(2a)^2<1>/\,(N_a\mathcal{N})\,(x,t)$','Interpreter','latex'); 
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
ylabel('$\overline{r_1}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,myy1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{r_2}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mVx1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{\xi_1}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mVy1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{\xi_2}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mxVx1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{r_1\xi_1}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mxVy1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{r_1\xi_2}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,myVx1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{r_2\xi_1}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,myVy1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{r_2\xi_2}(x,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

