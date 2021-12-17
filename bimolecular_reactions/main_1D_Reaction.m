%% CG averages over BGRW simulations of bimolecular reactions; computation of D & u; 
clear all; close all;
tic
%% Grid Initialization
it = 400;
I = it+1;
x1=0; x2=1;
dx = (x2-x1)/(I-1);
x = (x1:dx:x2);
T = 1;
tau=T/10;
a=(x2-x1)/20; ia=round(a/dx);
mt=(tau:4*tau:T);
mx=(3/2*a:0.1*a:x2-3/2*a); imx=round(mx/dx);
fprintf('a = %.2f tau = %.2f\n',a,tau);
%% Parameters
D1=0.01; % BGRW
% q=zeros(1,I);
q=1*ones(1,I);
Dfactor=2; % BGRW
dt=Dfactor*(2*D1/dx^2); dt=1./dt;
%% Initial Conditions
N=10^24; % = 1 mol; NI = nr. of particles / lattice point
NI1=round(N/I); NI2=round(N/I); % 1 mol uniformly distributed on I points
c10=NI1*ones(1,I); c20=NI2*ones(1,I);
c1BC=c10; c2BC=c20;
c1p=zeros(length(mt),I); c2p=zeros(length(mt),I);
mc1p=zeros(length(mt),length(mx)); mc2p=zeros(length(mt),length(mx));
mc1=zeros(length(mt),length(mx)); mc2=zeros(length(mt),length(mx));
mc1v=zeros(length(mt),length(mx)); mc2v=zeros(length(mt),length(mx));
mxx1=zeros(length(mt),length(mx)); mvv1=zeros(length(mt),length(mx)); mxv1=zeros(length(mt),length(mx));
mxx2=zeros(length(mt),length(mx)); mvv2=zeros(length(mt),length(mx)); mxv2=zeros(length(mt),length(mx));
%% Solution
t=0; 
while t<=T
    t=t+dt;
    %% Transport step
    nspec=2; %  molecular_species
    for i=1:nspec
        switch i
            case 1
                    [c1t,dmxx1,dmvv1,dmxv1]=CG_BGRW_1D(c10,c1BC,I,dx,dt,q,D1,mx,a); % BGRW transport solver
            case 2
                    [c2t,dmxx2,dmvv2,dmxv2]=CG_BGRW_1D(c20,c2BC,I,dx,dt,q,D1,mx,a); % BGRW transport solver
        end
    end
    %% reaction
    [c1,c2,dmc1,dmc2,dmc1p,dmc2p,dmc1v,dmc2v]=CG_reaction(c1t,c2t,dt,I,mx,imx,a,ia);
    c10=c1; c20=c2;
    %% CGST Averages
    for k=1:length(mt)
        if t>=(mt(k)-tau) && t<=(mt(k)+tau)
            if abs(t-mt(k))<dt/2
                c1p(k,:)=c1; c2p(k,:)=c2;
                mc1v(k,:)=dmc1v; mc2v(k,:)=dmc2v;
            end
            mc1p(k,:)=mc1p(k,:)+dmc1p*dt/(2*tau); mc2p(k,:)=mc2p(k,:)+dmc2p*dt/(2*tau);
            mc1(k,:)=mc1(k,:)+dmc1*dt/(2*tau); mc2(k,:)=mc2(k,:)+dmc2*dt/(2*tau);
            mxx1(k,:)=mxx1(k,:)+dmxx1*dt/(2*tau); mxx2(k,:)=mxx2(k,:)+dmxx2*dt/(2*tau);
            mvv1(k,:)=mvv1(k,:)+dmvv1*dt/(2*tau); mvv2(k,:)=mvv2(k,:)+dmvv2*dt/(2*tau);
            mxv1(k,:)=mxv1(k,:)+dmxv1*dt/(2*tau); mxv2(k,:)=mxv2(k,:)+dmxv2*dt/(2*tau);
            
            if round(t/dt)==round(mt(k)/dt)
                str=['(cg) t=',num2str(t)];
                strvect(k,1:length(str))=str;
            end
        end
    end
    tE=t;
end
%% Results
plot_CG2(strvect,x,mx,c1p,mc1,mc1p,mc1v,mxx1,mvv1,mxv1,c2p,mc2,mc2p,mc2v,D1,N);
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
fprintf('The total time is : %0.2e \n',tE) ;
fprintf('grid_Peclet = %0.4e \n',dx*mean(abs(q))/D1);
Dvect=-mxv1./mc1;
cgD=mean(mean(Dvect));
fprintf('cgD = %f \n',cgD);
sdcgD=std(mean(Dvect));
fprintf('SDcgD = %0.2e \n',sdcgD);
Dmax=cgD+sdcgD;
Dmin=cgD-sdcgD;
fprintf('cgDmax = %f cgDmin = %f \n',Dmax,Dmin);
Uvect=mvv1./mc1;
cgU=mean(mean(Uvect));
fprintf('cgU = %f \n',cgU);
sdcgU=std(mean(Uvect));
fprintf('SDcgU = %0.2e \n',sdcgU);
Umax=cgU+sdcgU;
Umin=cgU-sdcgU;
fprintf('cgUmax = %f cgUmin = %f \n',Umax,Umin);
toc
