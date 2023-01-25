%% Coupled degenerated flow & transport with Monod reactions
%% Deterministic-GRW flow solver; random-BGRW transport solver

clear all; 
close all
tic
global state;
initstate; state;
%%   Grid Initialization
% I=161; J=241; % for final results
I=41; J=61; % for tests
x1=0; x2=2;
y1=0; y2=3;
dx =(x2-x1)/(I-1);
x = x1:dx:x2;
dy=(y2-y1)/(J-1);
y=y1:dy:y2;
% T= 132;  % for final results
T= 1.1; % for tests
% tau=12; % for final results % 
tau=0.1; % for tests
a=(y2-y1)/30; 
ia=round(a/dx);
mt=[(T-tau)/3 2*(T-tau)/3 T-tau]; 
% imx=(I-1)/2; % centered sampling vertical profile
imx=7*(I-1)/8; % decentred
my=(3/2*a:0.5*a:y2-3/2*a); jmy=round(my/dx);
fprintf('a = %.2f tau = %.2f\n',a,tau);

%%   Parameters - loam soil
disp('silt loam')
Ksat = 4.96*10^-2;
theta_res=0.131;
theta_sat=0.396;
alpha=0.423;
nGM=2.06;
past=T/3;
t1=T/3;
maxr=0.8;
D1=0.001;
D2=D1;
Tolerance = 1e-5;
S= 50000;
Lp=10; Lc=1;
%% Initialization of random hydraulic conductivity K
Nmod = 100;
ZC1 = 0.1;
ZC2 = 0.001;
varK = 0.5;
[X,Y] = meshgrid(x,y);
[wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
C1 = wavenum(:,1);
C2 = wavenum(:,2);
Ks= K_r(X(:)',Y(:)',Nmod,Ksat,varK,C1,C2,phi);
Ks=reshape(Ks,J,I);N=10^24;
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p)  theta_GM(theta_res,theta_sat,p,alpha,nGM,ng);
K = @(tht) Ks.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2;
%% Initial Conditions
% IC - presure
p0=-((y'-y1)/(y2-y1))*3+0*x;
p = p0; pBC=p0;
% IC - concentrations
N=10^24;
i1=(I-1)/4+1; i2=3*(I-1)/4; 
dJ=(J-1)/10;
dw=i2-i1; % I;
NI1=round(N/(dw*dJ)); NI2=round(N/(I*J-dw*dJ));
c10 = zeros(J,I);
c20 = ones(J,I);
c10(J-dJ:J,i1:i2)=1;
c20(J-dJ:J,i1:i2)=0;
c10=c10*NI1; c20=c20*NI2;
cBC1=c10; cBC2=c20;
c1=c10; c2=c20;
%% figures Initial cnditions:
% figure;
% mesh(x,y,p);
% xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
% zlabel('$\psi(x,z,t=0)$','Interpreter','latex'); view(115,15);
% figure;
% mesh(x,y,c1/N);
% xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
% zlabel('$c_1(x,z,t=0)$','Interpreter','latex'); view(115,15);
% figure;
% mesh(x,y,c2/N);
% xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
% zlabel('$c_2(x,z,t=0)$','Interpreter','latex'); view(115,15);
% 
mc1=zeros(length(mt),length(my)); mc2=zeros(length(mt),length(my));
mc1v=zeros(length(mt),length(my)); mc2v=zeros(length(mt),length(my));
%% Solution
tgraf=0; t=0; kt=1;
Vx=zeros(J,I); Vy=zeros(J,I); pp=zeros(J,I);
tht = theta(p); tht0=tht;
thtc10=tht0.*c10; thtc20=tht0.*c20;
pa=p; c1a=c10; c2a=c20; 
Dfactor=1.2;
dtc=Dfactor*(2*D1/Lc/dx^2+2*D2/Lc/dy^2); dtc=1/dtc;
convf=zeros(3,S); convc1=zeros(3,S); convc2=zeros(3,S);
while t<=T
    eps=zeros(1,S); epsc1=zeros(1,S); epsc2=zeros(1,S);
    DK=K(tht);
    Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
    Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
    D=DK(2:J-1,2:I-1);
    dtp=4*(max(max(Dx))/(Lp*dx^2)+max(max(Dy))/(Lp*dy^2)); dtp=maxr*1/dtp;
    dt=min(dtc,dtp);
    t=t+dt;
    %% Flow step
    for s=1:S
        DK=K(tht);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        rx=dt*Dx/dx^2/Lp; ry=dt*Dy/dy^2/Lp;
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
            +ry(1:J-2,:).*p(1:J-2,2:I-1) +ry(2:J-1,:).*p(3:J,2:I-1);
        %% Boundary conditions - flow
        %%%% BCY Left/Right
        pp(:,1)=pp(:,2); % no flux left
        pp(:,I)=pp(:,I-1); % no flux right
        %%%% BCX Bottom/Top
        if t<=t1 % linear p(t) on \Gamma_1
            pp(J,:)=pBC(J,:)+2.2*t/t1;
        else
            pp(J,:)=0.2; % constant \psi(t) on \Gamma_1 % 
        end
%         pp(J,:)=pBC(J,:); % !!!! 
        %% Source term - flow
        dtht=(tht0-tht)/Lp;
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1);
        pp(2:J-1,2:I-1)=pp(2:J-1,2:I-1)+f;
        p=pp;
        tht = theta(p);
        p0=p;
        %% Convergence criterion - flow
        tol_eps=dx*norm(p-pa)+norm(p-pa)/norm(p);
        if kt*past>=t && kt*past<t+dt && t<=T
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end
        pa=p;
    end
    for s=1:S
        %% Transport step
        [Vx,Vy] = velocity(I,J,dx,dy,p0,D);
        thtc1=tht.*c1; thtc2=tht.*c2;
        nspec=2;
        for i=1:nspec
            switch i
                case 1
                    dthtc=(thtc10-thtc1)/Lc;
                    [c1t]=CG_BGRW_2D_L(c10,cBC1,dthtc,I,J,dx,dy,dt,Vx,Vy,D1,D2,Lc,dJ,i1,i2); % transport solver
                case 2
                    dthtc=(thtc20-thtc2)/Lc;
                    [c2t]=CG_BGRW_2D_L(c20,cBC2,dthtc,I,J,dx,dy,dt,Vx,Vy,D1,D2,Lc,dJ,i1,i2); % transport solver
            end
        end
        %% reaction
        [c1,c2,dmc1,dmc2,dmc1v,dmc2v]=CG_reaction_Monod2(c1t,c2t,dt,tht,Lc,I,J,my,imx,jmy,a,ia);
        c10=c1; c20=c2;
        %% Convergence criterion - transport
        tol_epsc1=norm(c1/N-c1a);
        tol_epsc2=norm(c2/N-c2a);
        if kt*past>=t && kt*past<t+dt && t<=T
            epsc1(s)=tol_epsc1;
            epsc2(s)=tol_epsc2;
        end
        if max(tol_epsc1,tol_epsc2) <= Tolerance
            break
        end
        c1a=c1/N; c2a=c2/N;
    end
    if  kt*past>=t && kt*past<t+dt && t<=T
        tgraf=tgraf+1;
        rndt=kt*past;
        fprintf('t = %0.2e \n',rndt) ;
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        convf(kt,:)=eps;
        convc1(kt,:)=epsc1;
        convc2(kt,:)=epsc2;
        kt=kt+1;
    end
    tht0=tht; thtc10=tht.*c1; thtc20=tht.*c2;
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
if  kt==4
    save('convf');
    save('convc1');
    save('convc2');
end
%% Results
kt_plot=kt-1;
plot_conv_Monod(kt_plot,S,strvect);
subplot_surf_c_2(x,y,c1,c2,N);
subplot_surf_ptht(x,y,p,tht,Vx,Vy)
plot_CG2(strvectm,my,mc1,mc1v,mc2,mc2v,N);

fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
%% mean and max Peclet numbers
V=mean(mean(sqrt(Vx.^2+Vy.^2)));
Pe=V*dx/D1;
maxV=max(max(sqrt(Vx.^2+Vy.^2)));
maxPe=maxV*dx/D1;
fprintf('meanV = %0.4e meanPe = %0.4e \n',V,Pe);
fprintf('maxV  = %0.4e maxPe  = %0.4e \n',maxV,maxPe);
toc
% Elapsed time is 222.771180 seconds = 3.7129 minutes for the test run
% = 5 hours for final results
