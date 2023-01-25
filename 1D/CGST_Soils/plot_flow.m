function plot_flow(k_flow,x,strvectm)

dk=10;
NameArray = {'Marker'}; ValueArray = {'o','+','*'}';

load pp
figure; hold all;
for k=1:k_flow
P(k)=plot(pp(k,1:dk:end),x(1:dk:end),'MarkerSize',4);
end
set(P,NameArray,ValueArray,LineWidth = 1); box on;
xlabel('$\psi$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
legend(strvectm,'Location','best'); legend('boxoff'); 

load thtp
figure; hold all;
for k=1:k_flow
P(k)=plot(thtp(k,1:dk:end),x(1:dk:end),'MarkerSize',4); box on;
end
set(P,NameArray,ValueArray,LineWidth = 1);
xlabel('$\theta$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
legend(strvectm,'Location','best'); legend('boxoff'); 

load qp
figure; hold all;
for k=1:k_flow
P(k)=plot(qp(k,1:dk:end),x(1:dk:end),'MarkerSize',4); box on;
end
set(P,NameArray,ValueArray,LineWidth = 1);
xlabel('$q$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
legend(strvectm,'Location','best'); legend('boxoff'); 
