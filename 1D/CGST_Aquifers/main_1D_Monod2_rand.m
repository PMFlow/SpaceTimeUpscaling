%% CG averages over BGRW simulations for bimolecular (Monod) reactions; stochastic aerages

clear all; close all;
tic
initstate; state;
%% Grid Initialization
irand=0; % 0=determinist; 1=random
R=1; % 100; % number of realizations
it = 400;
I = it+1;
x1=0; x2=1;
dx = (x2-x1)/(I-1);
x = (x1:dx:x2);
a=(x2-x1)/30; ia=round(a/dx);
T = 1.1; 
tau=0.1;
mt=[0.1 0.5 1];
mx=(1*a:0.1*a:x2-1*a); imx=round(mx/dx);
fprintf('a = %.2f tau = %.2f\n',a,tau);
%% Parameters
D1=0.01; 
Dfactor=2;
dt=Dfactor*(2*D1/dx^2); dt=1./dt;
c1p=zeros(length(mt),I); c2p=zeros(length(mt),I);
mc1p=zeros(length(mt),length(mx)); mc2p=zeros(length(mt),length(mx));
mc1=zeros(length(mt),length(mx)); mc2=zeros(length(mt),length(mx));
mc1v=zeros(length(mt),length(mx)); mc2v=zeros(length(mt),length(mx));

for r=1:R
% q=zeros(1,I);
q=1*ones(1,I);
norm_q=norm(q);
    if irand == 0 && R==1
        dq=zeros(1,I);
    else
        Nmod = 100;
        ZC1 = 0.01;
        ZC2 = ZC1;
        varq = 0.1;
        q0=1;
        [wavenum, phi, amplitude] = V_Kraichnan_Gauss_param(Nmod,varq,ZC1,ZC2,q0,0);
        for iq=1:I
            [dq(iq),vr] = V_Kraichnan_Gauss_func(x(iq),0.5,wavenum, phi, amplitude);    
        end
    end
    q=q+dq;
%% Initial Conditions
    N=10^24; % = 1 mole; 
% %   NI = nr. of particles / lattice point
% %    1 mole uniformly distributed on subdomains:
    dI=(I-1)/10;
    NI1=round(N/dI); NI2=round(N/(I-dI)); 
    c10=zeros(1,I); c20=NI2*ones(1,I);
    c10(1:dI)=NI1; c20(1:dI)=0;
    c1BC=c10; c2BC=c20;
    %% Solution
    t=0;
    while t<=T
        t=t+dt;
        %% Transport step
        nspec=2; %  molecular_species
        for i=1:nspec
            switch i
                case 1
                    [c1t]=CG_BGRW_1D(c10,c1BC,I,dI,dx,dt,q,D1); 
                case 2
                    [c2t]=CG_BGRW_1D(c20,c2BC,I,dI,dx,dt,q,D1); 
            end
        end
        %% reaction
        [c1,c2,dmc1,dmc2,dmc1p,dmc2p,dmc1v,dmc2v]=CG_reaction_Monod2(c1t,c2t,dt,I,mx,imx,a,ia,NI1,NI2);
        c10=c1; c20=c2;
        %% CGST Averages
        for k=1:length(mt)
            if t>=(mt(k)-tau) && t<=(mt(k)+tau)
                if abs(t-mt(k))<dt/2
                    c1p(k,:)=c1p(k,:)+c1; c2p(k,:)=c2p(k,:)+c2;
                    mc1v(k,:)=mc1v(k,:)+dmc1v; mc2v(k,:)=mc2v(k,:)+dmc2v;
                end
                mc1p(k,:)=mc1p(k,:)+dmc1p*dt/(2*tau); mc2p(k,:)=mc2p(k,:)+dmc2p*dt/(2*tau);
                mc1(k,:)=mc1(k,:)+dmc1*dt/(2*tau);    mc2(k,:)=mc2(k,:)+dmc2*dt/(2*tau);
                if round(t/dt)==round(mt(k)/dt)
                    str=['t=',num2str(t)];
                    strvect(k,1:length(str))=str;
                end
            end
        end
        tE=t;
    end
end
%% Results
plot_CG_Monod2(strvect,x,mx,c1p/R,mc1/R,mc1p/R,mc1v/R,c2p/R,mc2/R,mc2p/R,mc2v/R,N,R);
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
fprintf('The total time is : %0.2e \n',tE) ;
fprintf('grid_Peclet = %0.4e \n',dx*mean(abs(q))/D1);
toc

