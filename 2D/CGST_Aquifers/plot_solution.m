%% plot Kraichnan velocities and CG concentrations
function plot_solution(x,y,c1,c2,Vx,Vy)

N=10^24;

figure
mesh(x,y,Vy)
xlabel('$z$','Interpreter','latex'); ylabel('$x$','Interpreter','latex'); 
zlabel('$q_x(x,z)$','Interpreter','latex'); view(115,15);

figure
mesh(x,y,Vx)
xlabel('$z$','Interpreter','latex'); ylabel('$x$','Interpreter','latex'); 
zlabel('$q_z(x,z)$','Interpreter','latex'); view(115,15);

figure;
mesh(x,y,c1/N); 
xlabel('$z$','Interpreter','latex'); ylabel('$x$','Interpreter','latex');
zlabel('$c_1(x,z,t)$','Interpreter','latex'); view(115,15); 
grid on ;

figure;
mesh(x,y,c2/N); 
xlabel('$z$','Interpreter','latex'); ylabel('$x$','Interpreter','latex');
zlabel('$c_2(x,z,t)$','Interpreter','latex'); view(115,15); 
grid on ;

