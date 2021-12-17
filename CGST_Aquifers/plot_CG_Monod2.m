function plot_CG_Monod2(strvect,xx,x,c1,mc1,mc1p,mc1v,c2,mc2,mc2p,mc2v,N,R)

NameArray = {'Marker'}; ValueArray = {'o','+','*'}'; 

figure
P=plot(x(:,1:4:end),mc1p(:,1:4:end)/N,'MarkerSize',4); hold all; plot(xx,c1/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
if R>1
    ylabel('$M[2a<1_1>]/\,(N_a\mathcal{N})$','Interpreter','latex');
else
    ylabel('$2a<1_1>/\,(N_a\mathcal{N})$','Interpreter','latex'); 
end
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x(:,1:4:end),mc2p(:,1:4:end)/N,'MarkerSize',4); hold all; plot(xx,c2/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
if R>1
    ylabel('$M[2a<1_2>]/\,(N_a\mathcal{N})$','Interpreter','latex');
else
    ylabel('$2a<1_2>/\,(N_a\mathcal{N})$','Interpreter','latex');
end
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x(:,1:4:end),mc1(:,1:4:end)/N,'MarkerSize',4); hold all; plot(x,mc1v/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
if R>1
    ylabel('$M[<1_1>]/\,\mathcal{N}$','Interpreter','latex');
else
    ylabel('$<1_1>/\,\mathcal{N}$','Interpreter','latex');
end
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(x(:,1:4:end),mc2(:,1:4:end)/N,'MarkerSize',4); hold all; plot(x,mc2v/N,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1); 
xlabel('$x$','Interpreter','latex');
if R>1
    ylabel('$M[<1_2>]/\,\mathcal{N}$','Interpreter','latex');
else
    ylabel('$<1_2>/\,\mathcal{N}$','Interpreter','latex');
end
legend(strvect,'Location','best'); legend('boxoff');

