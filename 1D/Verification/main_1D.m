%% CG averages over BGRW/GRW simulations; computation of D & u

clear all; close all;
tic
%% Grid Initialization
I = 1001; % 2001; % 1601; % 801; % 401; % 201; % 
x1=0; x2=1; 
dx = (x2-x1)/(I-1);
x = (x1:dx:x2);
T = 1;
tau=T/10;
a=(x2-x1)/20; ia=round(a/dx);
mt=(tau:4*tau:T);
mx=(3/2*a:2.1*a:x2-3/2*a); imx=round(mx/dx);
fprintf('a = %.2f tau = %.2f\n',a,tau);
%% Parameters
igrw=1; % 1=BGRW (biased); 0=GRW (unbiased)
D1=0.0001; % GRW
q=1*ones(1,I); % constant velocity: 1* & 0*
norm_q=norm(q);
if igrw==1
    Dfactor=2; % BGRW
    dt=Dfactor*(2*D1/dx^2); dt=1./dt;
else
    d=1; stepU=1; U_MEAN=1; % GRW
        dt=stepU*dx/U_MEAN; 
end
%% Initial Conditions
N=10^24; % = 1 mol; NI = nr. of particles / lattice point
NI=round(N/I); % 1 mol uniformly distributed on I points: --->
% c_lattice point = 1/I mols & <1>/N=1/(2*a)*(mi/I);
% mi = nr. of points in [x-a, x+a]
c10=NI*ones(1,I); c1BC=c10;
mc1p=zeros(length(mt),length(mx)); mc1=zeros(length(mt),length(mx)); mc1v=zeros(length(mt),length(mx));
mxx1=zeros(length(mt),length(mx)); mvv1=zeros(length(mt),length(mx)); mxv1=zeros(length(mt),length(mx));
%% Solution
t=0;
while t<=T
    t=t+dt;
    %% Transport step
    if igrw==1
        [c1t,dmxx1,dmvv1,dmxv1]=CG_BGRW_1D(c10,I,dx,dt,q,norm_q,D1,mx,imx,a,ia); % BGRW transport solver
    else
        [c1t,dmxx1,dmvv1,dmxv1]=CG_GRW_1D(c10,I,dx,dt,q,D1,mx,imx,a,ia,d,stepU); % GRW transport solver
    end
    c1=c1t; c10=c1;   
    %% CGST Averages
    for k=1:length(mt)
        if t>(mt(k)-tau) && t<(mt(k)+tau)
            dmc1=zeros(1,length(mx));
            dmc1p=zeros(1,length(mx)); 
            dmc1v=zeros(1,length(mx));
            for m=1:length(mx)
                mi=0;
                for i=1:I
                    if abs(i-imx(m))<ia
                        dmc1(m)=dmc1(m)+c1(i)/(2*a);
                        dmc1v(m)=dmc1(m);
                        mi=mi+1;
                    end
                    dmc1p(m)=dmc1(m)*(2*a)/mi;
                end
%                 fprintf('mi= %d',mi);
            end
            mc1v(k,:)=dmc1v;
            mc1p(k,:)=mc1p(k,:)+dmc1p*dt/(2*tau);
            mc1(k,:)=mc1(k,:)+dmc1*dt/(2*tau);
            mxx1(k,:)=mxx1(k,:)+dmxx1*dt/(2*tau);
            mvv1(k,:)=mvv1(k,:)+dmvv1*dt/(2*tau);
            mxv1(k,:)=mxv1(k,:)+dmxv1*dt/(2*tau);

            if round(t/dt)==round(mt(k)/dt)
                str=['t=',num2str(t)];
                strvect(k,1:length(str))=str;
                ntot=sum(c1);
            end
        end
    end
    tE=t;
end
%% Results
plot_CG(strvect,x,mx,c1,mc1,mc1p,mc1v,mxx1,mvv1,mxv1,D1,N);
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
fprintf('The total time is : %0.2e \n',tE) ;
fprintf('grid_Peclet = %0.4e \n',dx*mean(abs(q))/D1);
Dvect=-mxv1./mc1;
cgD=mean(mean(Dvect));
sdcgD=std(mean(Dvect));
Dmax=cgD+sdcgD;
Dmin=cgD-sdcgD;
Uvect=mvv1./mc1;
cgU=mean(mean(Uvect));
sdcgU=std(mean(Uvect));
Umax=cgU+sdcgU;
Umin=cgU-sdcgU;
fprintf('cgD = %0.2e +/- %0.2e\n',cgD,sdcgD);
fprintf('cgU = %0.2f +/- %0.2e\n',cgU,sdcgU);
toc
