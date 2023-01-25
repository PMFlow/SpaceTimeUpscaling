function plot_CG2(strvect,xx,x,c1,mc1,mc1p,mc1v,mxx1,mvv1,mxv1,c2,mc2,mc2p,mc2v,D,N)

NameArray = {'Marker'}; ValueArray = {'o','+','*'}'; 
figure
P=plot(x(:,1:4:end),mc1p(:,1:4:end)/N,'MarkerSize',4); hold all; plot(xx,c1/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
ylabel('$2a<1_1>/\,(N_a\mathcal{N})$','Interpreter','latex'); 
plot(xx(:,1:20:end),(c1(:,1:20:end)+c2(:,1:20:end))/2/N,'ks-','MarkerSize',4,LineWidth = 0.75);
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x(:,1:4:end),mc2p(:,1:4:end)/N,'MarkerSize',4); hold all; plot(xx,c2/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
ylabel('$2a<1_2>/\,(N_a\mathcal{N})$','Interpreter','latex');
plot(xx(:,1:20:end),(c1(:,1:20:end)+c2(:,1:20:end))/2/N,'ks-','MarkerSize',4,LineWidth = 0.75);
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x(:,1:4:end),mc1(:,1:4:end)/N,'MarkerSize',4); hold all; plot(x,mc1v/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
ylabel('$<1_1>/\,\mathcal{N}$','Interpreter','latex');
plot(x(:,1:10:end),(mc1(:,1:10:end)+mc2(:,1:10:end))/2/N,'k-s','MarkerSize',4,LineWidth = 0.75); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x(:,1:4:end),mc2(:,1:4:end)/N,'MarkerSize',4); hold all; plot(x,mc2v/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
ylabel('$<1_2>/\,\mathcal{N}$','Interpreter','latex');
plot(x(:,1:10:end),(mc1(:,1:10:end)+mc2(:,1:10:end))/2/N,'k-s','MarkerSize',4,LineWidth = 0.75); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mxx1./mc1); 
set(P,NameArray,ValueArray,LineWidth = 1.5); 
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{x}_1(z,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mvv1./mc1); 
set(P,NameArray,ValueArray,LineWidth = 1.5); 
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{\xi}_1(z,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,mxv1./mc1);
set(P,NameArray,ValueArray,LineWidth = 1.5); 
xlabel('$x$','Interpreter','latex');
ylabel('$\overline{x\xi}_1(z,t)$','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x,(-mxv1./mc1/D-1)*100);
set(P,NameArray,ValueArray,LineWidth = 1.5);  
xlabel('$x$','Interpreter','latex');
ylabel('$(cg\_D-D)/D$ [\%]','Interpreter','latex'); 
legend(strvect,'Location','best'); legend('boxoff');
 