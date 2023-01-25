%% 2-dimensional CG averages over BGRW/GRW simulations; computation of Dij, u & v

clear all; close all
tic
%%   Grid Initialization
I=201; % 401; % 801; % 1601; 
J=I;
x1=0; x2=1;
y1=0; y2=1;

dx = (x2-x1)/(I-1);
x = x1:dx:x2;
dy=(y2-y1)/(J-1);
y=y1:dy:y2;
T = 1;
tau=T/10;
a=(x2-x1)/20; ia=round(a/dx);
mt=(tau:4*tau:T);
mx=(3/2*a:2.1*a:x2-3/2*a); imx=round(mx/dx);
my=(J-1)/2*dy; jmy=(J-1)/2;
fprintf('a = %.2f tau = %.2f\n',a,tau);
%% Parameters
igrw=1; % 1=BGRW (biased); 0=GRW (unbiased)
D1=0.0001; D2=D1; 
Vx=ones(J,I); % constant velocity 
Vy=zeros(J,I); 
normVx=norm(Vx); normVy=norm(Vy);
if igrw==1
    Dfactor=2; % BGRW
    dt=Dfactor*(2*D1/dx^2+2*D2/dy^2); dt=1/dt;
else
    d=1; stepU=1; U_MEAN=1; % GRW
    dt=stepU*dx/U_MEAN; 
end
%% Initial Conditions
N=10^24; % = 1 mole; NI = nr. of particles / lattice point
NI=round(N/(I*J)); % 1 mole uniformly distributed on IxJ points: ---> 
% c_lattice point = (1/N)*NI = 1/(I*J) moles & <1>/N=1/((2*a)^2)*(mi/(I*J)) = space averge; 
% mi = nr. of points in [x-a, x+a]^2
c10=NI*ones(J,I); c1BC=c10;
mc1p=zeros(length(mt),length(mx)); mc1=zeros(length(mt),length(mx)); mc1v=zeros(length(mt),length(mx));
mxx1=zeros(length(mt),length(mx)); myy1=zeros(length(mt),length(mx));
mVx1=zeros(length(mt),length(mx)); mVy1=zeros(length(mt),length(mx));
mxVx1=zeros(length(mt),length(mx)); mxVy1=zeros(length(mt),length(mx));
myVx1=zeros(length(mt),length(mx)); myVy1=zeros(length(mt),length(mx));
%% Solution
t=0; 
while t<=T
    t=t+dt;
    %% Transport step
    if igrw==1
        [c1t,dmxx1,dmyy1,dmVx1,dmVy1,dmxVx1,dmxVy1,dmyVx1,dmyVy1]=CG_BGRW_2D(c10,I,J,dx,dy,dt,Vx,Vy,normVx,normVy,D1,D2,mx,a,ia,imx,jmy); % BGRW transport solver
    else
        [c1t,dmxx1,dmyy1,dmVx1,dmVy1,dmxVx1,dmxVy1,dmyVx1,dmyVy1]=CG_GRW_2D(c10,I,J,dx,dy,dt,Vx,Vy,D1,D2,mx,a,ia,imx,jmy,d,stepU); % GRW transport solver
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
                for j=1:J
                    for i=1:I
                        if abs(i-imx(m))<ia && abs(j-jmy)<ia
                            dmc1(m)=dmc1(m)+c1(j,i)/(2*a)^2;
                            dmc1v(m)=dmc1(m);
                            mi=mi+1;
                        end
                        dmc1p(m)=dmc1(m)*(2*a)^2/mi;
                    end
                end
            end
            mc1v(k,:)=dmc1v;
            mc1p(k,:)=mc1p(k,:)+dmc1p*dt/(2*tau);
            mc1(k,:)=mc1(k,:)+dmc1*dt/(2*tau);
            mxx1(k,:)=mxx1(k,:)+dmxx1*dt/(2*tau); myy1(k,:)=myy1(k,:)+dmyy1*dt/(2*tau);
            mVx1(k,:)=mVx1(k,:)+dmVx1*dt/(2*tau); mVy1(k,:)=mVy1(k,:)+dmVy1*dt/(2*tau);
            mxVx1(k,:)=mxVx1(k,:)+dmxVx1*dt/(2*tau); mxVy1(k,:)=mxVy1(k,:)+dmxVy1*dt/(2*tau);
            myVx1(k,:)=myVx1(k,:)+dmyVx1*dt/(2*tau); myVy1(k,:)=myVy1(k,:)+dmyVy1*dt/(2*tau);

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
plot_CG(strvect,x,mx,c1(jmy,:),mc1,mc1p,mc1v,mxx1,myy1,mVx1,mVy1,mxVx1,mxVy1,myVx1,myVy1,N);
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
fprintf('The total time is : %0.2e \n',tE) ;

q=sqrt(Vx.^2+Vy.^2);
fprintf('grid_Peclet = %0.4e \n',dx*mean(mean(q))/D1);

Dvect=-mxVx1./mc1;
cgD=mean(mean(Dvect));
sdcgD=std(mean(Dvect));
fprintf('cgD_{11} = %0.2e +/- %0.2e\n',cgD,sdcgD);

Dvect=-mxVy1./mc1;
cgD=mean(mean(Dvect));
sdcgD=std(mean(Dvect));
fprintf('cgD_{12} = %0.2e +/- %0.2e\n',cgD,sdcgD);

Dvect=-myVx1./mc1;
cgD=mean(mean(Dvect));
sdcgD=std(mean(Dvect));
fprintf('cgD_{21} = %0.2e +/- %0.2e\n',cgD,sdcgD);

Dvect=-myVy1./mc1;
cgD=mean(mean(Dvect));
sdcgD=std(mean(Dvect));
fprintf('cgD_{22} = %0.2e +/- %0.2e\n',cgD,sdcgD);

Uvect=mVx1./mc1;
cgU=mean(mean(Uvect));
sdcgU=std(mean(Uvect));
fprintf('cgV_1 = %0.2e +/- %0.2e\n',cgU,sdcgU);

Uvect=mVy1./mc1;
cgU=mean(mean(Uvect));
sdcgU=std(mean(Uvect));
fprintf('cgV_2 = %0.2e +/- %0.2e\n',cgU,sdcgU);
toc

%--------------------------------------------------------------------------

%% main_2D: BGRW, I=J=1601
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 6.25e-04 
% The time step is : 4.88e-04 
% The total time is : 1.00e+00 
% grid_Peclet = 6.2500e+00 
% cgD_{11} = 1.00e-04 +/- 7.47e-19
% cgD_{12} = 1.23e-19 +/- 3.02e-19
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 1.44e-20
% cgV_1 = 1.00e+00 +/- 2.36e-16
% cgV_2 = -9.55e-20 +/- 1.28e-35
% Elapsed time is 1781.273849 seconds = 29.6879 min.

%% main_2D: BGRW, I=J=801
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 1.25e-03 
% The time step is : 1.95e-03 
% The total time is : 1.00e+00 
% grid_Peclet = 1.2500e+01 
% cgD_{11} = 1.00e-04 +/- 3.02e-19
% cgD_{12} = 5.93e-20 +/- 2.36e-19
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 2.87e-20
% cgV_1 = 1.00e+00 +/- 0.00e+00
% cgV_2 = 1.14e-19 +/- 2.55e-35
% Elapsed time is 110.112055 seconds = 1.8352 min.

%% main_2D: BGRW, I=J=401
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 2.50e-03 
% The time step is : 7.81e-03 
% The total time is : 1.01e+00 
% grid_Peclet = 2.5000e+01 
% cgD_{11} = 1.00e-04 +/- 1.31e-19
% cgD_{12} = -1.06e-19 +/- 2.92e-19
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 0.00e+00
% cgV_1 = 1.00e+00 +/- 1.18e-16
% cgV_2 = -1.39e-19 +/- 0.00e+00
% Elapsed time is 7.740926 seconds.

%% main_2D: BGRW, I=J=201
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 5.00e-03 
% The time step is : 3.12e-02 
% The total time is : 1.03e+00 
% grid_Peclet = 5.0000e+01 
% cgD_{11} = 1.00e-04 +/- 1.59e-19
% cgD_{12} = -1.59e-20 +/- 8.13e-20
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 0.00e+00
% cgV_1 = 1.00e+00 +/- 0.00e+00
% cgV_2 = -2.71e-19 +/- 0.00e+00
% Elapsed time is 2.190359 seconds.

%--------------------------------------------------------------------------

%% main_2D: GRW, I=J=1601
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 6.25e-04 
% The time step is : 6.25e-04 
% The total time is : 1.00e+00 
% grid_Peclet = 6.2500e+00 
% cgD_{11} = 1.00e-04 +/- 9.66e-19
% cgD_{12} = 2.07e-19 +/- 3.56e-19
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 1.44e-20
% cgV_1 = 1.00e+00 +/- 0.00e+00
% cgV_2 = 0.00e+00 +/- 0.00e+00
% Elapsed time is 1233.367099 seconds = 20.5561 min.

%% main_2D: GRW, I=J=801
% main_2D
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 1.25e-03 
% The time step is : 1.25e-03 
% The total time is : 1.00e+00 
% grid_Peclet = 1.2500e+01 
% cgD_{11} = 1.00e-04 +/- 2.46e-19
% cgD_{12} = -1.09e-19 +/- 2.73e-19
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 1.44e-20
% cgV_1 = 1.00e+00 +/- 0.00e+00
% cgV_2 = 0.00e+00 +/- 0.00e+00
% Elapsed time is 146.157273 seconds = 2.4360 min.

%% main_2D: GRW, I=J=401
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 2.50e-03 
% The time step is : 2.50e-03 
% The total time is : 1.00e+00 
% grid_Peclet = 2.5000e+01 
% cgD_{11} = 1.00e-04 +/- 1.63e-19
% cgD_{12} = -5.21e-20 +/- 2.63e-19
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 2.87e-20
% cgV_1 = 1.00e+00 +/- 0.00e+00
% cgV_2 = 0.00e+00 +/- 0.00e+00
% Elapsed time is 18.516009 seconds.

%% main_2D: GRW, I=J=201
% a = 0.05 tau = 0.10
% Current plot held
% The space step is : 5.00e-03 
% The time step is : 5.00e-03 
% The total time is : 1.00e+00 
% grid_Peclet = 5.0000e+01 
% cgD_{11} = 1.00e-04 +/- 1.20e-19
% cgD_{12} = -3.57e-21 +/- 7.29e-20
% cgD_{21} = -0.00e+00 +/- 0.00e+00
% cgD_{22} = 1.00e-04 +/- 1.44e-20
% cgV_1 = 1.00e+00 +/- 0.00e+00
% cgV_2 = 0.00e+00 +/- 0.00e+00
% Elapsed time is 3.624446 seconds.

%--------------------------------------------------------------------------
