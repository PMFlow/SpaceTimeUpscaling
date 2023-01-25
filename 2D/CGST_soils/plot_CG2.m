%% CGST averages along veritcal sampling lines

function plot_CG2(strvect,x,mc1,mc1v,mc2,mc2v,N)

NameArray = {'Marker'}; ValueArray = {'o','+','*'}'; 

figure
P=plot(mc1(:,1:1:end)/N,x(:,1:1:end),'MarkerSize',4,LineWidth = 1); hold all; plot(mc1v/N,x,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1);  
xlabel('$<1_1>/\,\mathcal{N}$','Interpreter','latex'); 
ylabel('$z$','Interpreter','latex');
legend(strvect,'Location','best'); legend('boxoff');

figure
P=plot(mc2(:,1:1:end)/N,x(:,1:1:end),'MarkerSize',4,LineWidth = 1); hold all; plot(mc2v/N,x,'g-',LineWidth = 1); 
set(P,NameArray,ValueArray,LineWidth = 1);  
xlabel('$<1_2>/\,\mathcal{N}$','Interpreter','latex'); 
ylabel('$z$','Interpreter','latex');
legend(strvect,'Location','best'); legend('boxoff');
