%% CG averages over BGRW simulations of bimolecular (Monod) reactions in soils
clear all; close all;
tic
initstate; state;
%% Grid Initialization
it = 240; 
I = it+1;
x1=0; x2=3; 
dx = (x2-x1)/(I-1);
x = (x1:dx:x2);
T = 6; 
tau=T/10; a=(x2-x1)/60;
ia=round(a/dx);
mt=(tau:4*tau:T-tau);
mx=(3/2*a:1*a:x2-3/2*a); imx=round(mx/dx);
fprintf('a = %.2f tau = %.2f\n',a,tau);
past=T/3;
t1=T/3;
Tolerance = 1e-5;
S= 50000;
maxr=1;
%% Soil parameters
disp('slit loam')
Ksat = 4.96*10^-2;
theta_res=0.131;
theta_sat=0.396;
alpha=0.423;
nGM=2.06;
Lp=10; Lc=Lp;
%% Diffusion parameters
D1=0.001; 
Dfactor=2; % BGRW
dtc=Dfactor*(2*D1/Lc/dx^2); dtc=1./dtc; 
%% Initialization of random/deterministic hydraulic conductivity K
Nmod = 100;
ZC1 = 0.1;
ZC2 = 0.001;
varK = 0.5; 
[wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
C1 = wavenum(:,1);
C2 = wavenum(:,2);
Ks= K_r(x(:)',Nmod,Ksat,varK,C1,C2,phi);
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p)  theta_GM_flow(theta_res,theta_sat,p,alpha,nGM,ng);
K = @(tht) Ks.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2;
%%   Initial conditions - pressure
p0=1-((x-x1)/(x2-x1))*3;
p = p0; pBC=p0;
%% Initial conditions -comcetrations
dI=(I-1)/10;
N=10^24; % = 1 mol; NI = nr. of particles / lattice point
NI1=round(N/dI); NI2=round(N/(I-dI)); 
c10=zeros(1,I); c20=ones(1,I); 
c10(I-dI:I)=1; c20(I-dI:I)=0; 
c10=c10*NI1; c20=c20*NI2;
c1BC=c10; c2BC=c20;
c1=c10; c2=c20; 
pp=zeros(length(mt),I); thtp=zeros(length(mt),I); qp=zeros(length(mt),I); 
c1p=zeros(length(mt),I); c2p=zeros(length(mt),I);
mc1p=zeros(length(mt),length(mx)); mc2p=zeros(length(mt),length(mx));
mc1=zeros(length(mt),length(mx)); mc2=zeros(length(mt),length(mx));
mc1v=zeros(length(mt),length(mx)); mc2v=zeros(length(mt),length(mx));
%% Solution
convf=zeros(3,S); convc1=zeros(3,S); convc2=zeros(3,S);
tht=theta(p); tht0=tht;
thtc10=tht0.*c10; thtc20=tht0.*c20;
pa=p; c1a=c10; c2a=c20;
n0 = floor(N*p); n0E=round(pBC*N);
nn=zeros(I); n=n0;
restr=0; restsar1=0; restsarI=0; rest1=0; rest2=0; restf=0;
t=0; dt=dtc; tgraf=0; kt=1;
while t<=T
    eps=zeros(1,S); epsc1=zeros(1,S); epsc2=zeros(1,S);
    DK=K(tht);
    D=(DK(1:I-1)+DK(2:I))/2;
    Dq=DK(2:I-1);
    dtp=Lp*dx^2*maxr/max(D)/2; 
    dt=min(dtc,dtp); 
    t=t+dt;
%% Flow solution
    for s=1:S
        %% Pressure step
        DK=K(tht);
        D=(DK(1:I-1)+DK(2:I))/2;
        Dq=DK(2:I-1);
        r=dt*D/dx^2/Lp;
        rloc=[1-2*r(1),1-(r(1:I-2)+r(2:I-1)),1-2*r(I-1)];
        rapr=[0.5,r(1:I-2)./(r(1:I-2)+r(2:I-1))];
        restr=rloc.*n+restr; nn=floor(restr); restr=restr-nn;  nsar=n-nn;
        restsar1=r(1)*n(1)+restsar1; nsar1=floor(restsar1); restsar1=restsar1-nsar1;
        nn(2)=nn(2)+nsar1;
        rest1=r(1:I-2).*n(2:I-1)+rest1; nsarleft=floor(rest1); rest1=rest1-nsarleft;
        nn(1:I-2)=nn(1:I-2)+nsarleft;
        nn(3:I)=nn(3:I)+nsar(2:I-1)-nsarleft;
        restsarI=r(I-1)*n(I)+restsarI; nsarI=floor(restsarI); restsarI=restsarI-nsarI;
        nn(I-1)=nn(I-1)+nsarI;
        %% Boundary conditions - pressure
        nn(1)=n0E(1);
        if t<=t1
            nn(I)=n0E(I)+2.2*t/t1;
        else
            nn(I)=0.2;
        end
        %% Source term - pressure
        dtht=(tht0-tht)/Lp;
        f=diff(r)*dx+dtht(2:I-1);
        restf=N*f+restf; nf=floor(restf); restf=restf-nf;
        nn(2:I-1)=nn(2:I-1)+nf;
        n=nn; p=n/N;
        p0=p; tht=theta(p);
        %% Convergence criterion - pressure
        tol_eps=dx*norm(p-pa)+norm(p-pa)/norm(p);
        if kt*past>=t && kt*past<t+dt && t<=T
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end
        pa=p; 
    end
    %% Velocity step
    q(2:I-1)=-Dq.*((p0(3:I)-p0(1:I-2))/(2*dx)+1); % flow-interior of \Omega
    q(1)=-DK(1).*((p0(2)-p0(1))/dx+1); % flow left_BC approximated by finie-difference
    q(I)=-DK(I).*((p0(I)-p0(I-1))/dx+1); % flow right_BC approximated by finie-difference
    mean_q=mean(q)*ones(1,I);
    for s=1:S
    %% Transport step
    thtc1=tht.*c1; thtc2=tht.*c2;
    nspec=2; %  molecular_species
    for i=1:nspec
        switch i
            case 1
                dthtc=(thtc10-thtc1)/Lc;
                [c1t]=CG_BGRW_1D_L(c10,c1BC,dthtc,I,dI,dx,dt,q,D1,Lc); % BGRW transport solver
            case 2
                dthtc=(thtc20-thtc2)/Lc;
                [c2t]=CG_BGRW_1D_L(c20,c2BC,dthtc,I,dI,dx,dt,q,D1,Lc); % BGRW transport solver
        end
    end
    %% reaction
    [c1,c2,dmc1,dmc2,dmc1p,dmc2p,dmc1v,dmc2v]=CG_reaction_Monod2(c1t,c2t,dt,I,mx,imx,a,ia,NI1,NI2);
    c10=c1; c20=c2;
        %% Convergence criterion - transport
        tol_epsc1=norm(c1/NI1-c1a);
        tol_epsc2=norm(c2/NI2-c2a);
        if kt*past>=t && kt*past<t+dt && t<=T
            epsc1(s)=tol_epsc1;
            epsc2(s)=tol_epsc2;
        end
        if max(tol_epsc1,tol_epsc2) <= Tolerance
            break
        end
        c1a=c1/NI1; c2a=c2/NI2;
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
    tE=t;
    %% CGST Averages
    for k=1:length(mt)
        if t>=(mt(k)-tau) && t<=(mt(k)+tau)
            if abs(t-mt(k))<dt/2
                c1p(k,:)=c1; c2p(k,:)=c2;
                mc1v(k,:)=dmc1v; mc2v(k,:)=dmc2v;
            end
            mc1p(k,:)=mc1p(k,:)+dmc1p*dt/(2*tau); mc2p(k,:)=mc2p(k,:)+dmc2p*dt/(2*tau);
            mc1(k,:)=mc1(k,:)+dmc1*dt/(2*tau); mc2(k,:)=mc2(k,:)+dmc2*dt/(2*tau);
            if round(t/dt)==round(mt(k)/dt)
                strm=['t=',num2str(mt(k))]; 
                strvectm(k,1:length(strm))=strm;
                pp(k,:)=p;
                thtp(k,:)=tht;
                qp(k,:)=q;
            end
        end
    end
end
if  kt==4
    save('convf');
    save('convc1');
    save('convc2');
end
save('pp');
save('thtp');
save('qp');

%% Results
kt_plot=kt-1;
plot_conv(kt_plot,S,strvect);
k_flow=length(mt);
plot_flow(k_flow,x,strvectm)
plot_CG2(strvectm,x,mx,c1p,mc1,mc1p,mc1v,c2p,mc2,mc2p,mc2v,N);
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
fprintf('The total time is : %0.2e \n',tE) ;
fprintf('grid_Peclet = %0.4e \n',dx*mean(abs(q))/D1);

toc
