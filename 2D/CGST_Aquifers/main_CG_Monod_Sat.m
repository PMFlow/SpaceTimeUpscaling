%% Reactive transport in saturated aquifers
%% Kraichnan velocities; biased GRW transport solver; 
%% Monod reactions; electron acceptor and electron donor; constant biomass concetration
clear all; close all
tic
global state;
initstate; state;
irand=1; % 0=constant Ksat; 1=random Ksat;
%%   Grid Initialization
I=101; J=201;  
x1=0; x2=10; 
y1=0; y2=20; 

dx = (x2-x1)/(I-1);
x = x1:dx:x2;
dy=(y2-y1)/(J-1);
y=y1:dy:y2;
%%
T= 330; 
tau=T/11;
a=(y2-y1)/30; 
ia=round(a/dx);
mt=[(T-tau)/3 2*(T-tau)/3 T-tau]; 
imx=(I-1)/2; % sampling on centered vertical profile
my=(3/2*a:0.5*a:y2-3/2*a); jmy=round(my/dx);
fprintf('a = %.2f tau = %.2f\n',a,tau);

%%   Transport parameters
Ksat = 2; 
tht = 1; 
D1=0.01; D2=D1; 
stepU=2; 
p0 = 0.8 - 0.8*((y'-y1)/(y2-y1))+0*x; 
U_MEAN=-Ksat*((p0(J,2)-p0(1,2))/(y2-y1))/max(max(tht));
fprintf('MeanVx = %.2f\n',U_MEAN);
Dfactor=1.2;
dtc=Dfactor*(2*D1/dx^2/min(min(tht))+2*D2/dy^2/min(min(tht))); dtc=1./dtc; % BGRW
dt=dtc;

%% Initialization Karaichnan routine
NMOD=100; % number of random modes
K_MEAN=Ksat; % meand K-field
varK=0.5; % variance of lnK-field
ZC1=1.0; ZC2=1.0; % correlation lengths
lambda=0.0; % width of spatial filtering
[X,Y] = meshgrid(x,y);
[wavenum, phi, amplitude] = V_Kraichnan_Gauss_param(NMOD,varK,ZC1,ZC2,U_MEAN,lambda);
    Vx=zeros(J,I); Vy=zeros(J,I);   
    for j=1:J
        for i=1:I
            if irand==0
                ur=0; 
                vr=0;
            else
                yy=j*dy; xx=i*dx; 
                [ur,vr] = V_Kraichnan_Gauss_func(xx,yy,wavenum, phi, amplitude);
            end
            Vx(j,i)=ur; 
            Vy(j,i)=vr+U_MEAN; 
        end
    end
%%   Initial Conditions - concentrations   %% 
N=10^24;
c10 = zeros(J,I);
c20 = 0.1*ones(J,I); 
i0=round(I/2); j0=round(1/dx); dw=round(0.5/dx); 
c10(j0-dw:j0+dw,i0-dw:i0+dw)=2.5; 
c20(j0-dw:j0+dw,i0-dw:i0+dw)=0;
c10=c10*N; c20=c20*N;  
cBC1=c10; cBC2=c20; 
c1=c10; c2=c20; 
% figure;
% mesh(x,y,c1);
% xlabel('z','Interpreter','latex'); ylabel('x','Interpreter','latex');
% zlabel('$c_1(x,z,t=0)$','Interpreter','latex'); view(115,15);
% figure;
% mesh(x,y,c2);
% xlabel('z','Interpreter','latex'); ylabel('x','Interpreter','latex');
% zlabel('$c_2(x,z,t=0)$','Interpreter','latex'); view(115,15);

%% Solution
t=0; 
mc1=zeros(length(mt),length(my)); mc2=zeros(length(mt),length(my));
mc1v=zeros(length(mt),length(my)); mc2v=zeros(length(mt),length(my));
%% Transport step
while t<=T
    t=t+dt;
    nspec=2;
    for i=1:nspec
        switch i
            case 1
                [c1t]=BGRW_2D_React_Monod_Sat(c10,cBC1,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,dw); % transport solver
            case 2
                [c2t]=BGRW_2D_React_Monod_Sat(c20,cBC2,tht,I,J,dx,dy,i0,j0,dt,Vx,Vy,D1,D2,dw); % transport solver
        end
    end
    %% reaction
    [c1,c2,dmc1,dmc2,dmc1v,dmc2v]=CG_reaction_Monod2(c1t,c2t,dt,tht,I,J,my,imx,jmy,a,ia);
    c10=c1; c20=c2;
    %% CGST Averages
    for k=1:length(mt)
        if t>(mt(k)-tau) && t<(mt(k)+tau)
            mc1(k,:)=mc1(k,:)+dmc1*dt/(2*tau); mc2(k,:)=mc2(k,:)+dmc2*dt/(2*tau);
            if abs(t-mt(k))<dt/2
                mc1v(k,:)=dmc1v; mc2v(k,:)=dmc2v;
                strm=['t=',num2str(t)];
                strvectm(k,1:length(strm))=strm;
            end
        end
    end
    
end

%% Results
plot_solution(x,y,c1,c2,Vx,Vy)
plot_surf_c(y,x,c1',c2',N);
plot_CG2(strvectm,my,mc1,mc1v,mc2,mc2v,N);
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
%% mean and max Peclet numbers
V=mean(mean(sqrt(Vx.^2+Vy.^2)));
Pe=V*dx/D1;
maxV=max(max(sqrt(Vx.^2+Vy.^2)));
maxPe=maxV*dx/D1;
fprintf('mean|V| = %0.4e meanPe = %0.4e \n',V,Pe);
fprintf('max|V|  = %0.4e maxPe  = %0.4e \n',maxV,maxPe);
errors_aqv(mc1,mc1v,mc2,mc2v)
toc 
