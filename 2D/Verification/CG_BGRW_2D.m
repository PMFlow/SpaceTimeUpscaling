function [c,dmxx,dmyy,dmVx,dmVy,dmxVx,dmxVy,dmyVx,dmyVy]=CG_BGRW_2D(c0,I,J,dx,dy,dt,Vx,Vy,normVx,normVy,D1,D2,mx,a,ia,imx,jmy)

%% BGRW solution
v=Vy*dt/dx; u=Vx*dt/dy;
ru=2*D1*dt/dx^2*ones(J,I);
rv=2*D2*dt/dy^2*ones(J,I);
n=c0;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0; restjgrw=0;
dmxx=zeros(1,length(mx)); dmyy=zeros(1,length(mx));
dmVx=zeros(1,length(mx)); dmVy=zeros(1,length(mx));
dmxVx=zeros(1,length(mx)); dmxVy=zeros(1,length(mx));
dmyVx=zeros(1,length(mx)); dmyVy=zeros(1,length(mx));

for y=1:J
    for x=1:I
        if n(y,x) > 0
            rx=ru(y,x); ry=rv(y,x); r=rx+ry;
            restr=n(y,x)*(1-r)+restr; nsta=floor(restr);
            restr=restr-nsta; njump=n(y,x)-nsta;
            nn(y,x)=nn(y,x)+nsta;
            restjump=njump*ry/r+restjump;
            njumpy=floor(restjump); restjump=restjump-njumpy;
            njumpx=njump-njumpy;
            for m=1:length(mx)
                if abs(x-imx(m))<ia && abs(y-jmy)<ia
                    dmxx(m)=dmxx(m)+nsta*x*dx/(2*a)^2;
                    dmyy(m)=dmyy(m)+nsta*y*dy/(2*a)^2;
                end
            end
            if(njumpy)>0
                restjumpy=njumpy*0.5*(1-v(y,x)/ry)+restjumpy;
                nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                nj(2)=njumpy-nj(1);
                %
                restjgrw=njumpy/2+restjgrw; % ~unbiased GRW
                njgrw(1)=floor(restjgrw); restjgrw=restjgrw-njgrw(1);
                njgrw(2)=njumpy-njgrw(1);
                %
                if y==1
                    nn(2,x)=nn(2,x)+nj(2);
                elseif y==J
                    nn(J-1,x)=nn(J-1,x)+nj(1);
                else
                    for i=1:2
                        yd=y+(2*i-3);
                        nn(yd,x)=nn(yd,x)+nj(i);
                        for m=1:length(mx)
                            if abs(x-imx(m))<ia && abs(y-jmy)<=ia && abs(yd-jmy)<ia
                                dmyy(m)=dmyy(m)+nj(i)*yd*dy/(2*a)^2;
                                dmxx(m)=dmxx(m)+nj(i)*x*dx/(2*a)^2;
                                dmVy(m)=dmVy(m)+nj(i)*(yd-y)*dy/dt/(2*a)^2;
                                if normVy==0
                                    dmyVy(m)=dmyVy(m)+nj(i)*(y+yd)/2*dy*(yd-y)*dy/dt/(2*a)^2;
                                    dmxVy(m)=dmxVy(m)+nj(i)*x*dx*(yd-y)*dy/dt/(2*a)^2;
                                else
                                    dmyVy(m)=dmyVy(m)+njgrw(i)*(y+yd)/2*dy*(yd-y)*dy/dt/(2*a)^2;
                                    dmxVy(m)=dmxVy(m)+njgrw(i)*x*dx*(yd-y)*dy/dt/(2*a)^2;
                                end
                            end
                        end
                    end
                end
            end
            if(njumpx)>0
                restjumpx=njumpx*0.5*(1-u(y,x)/rx)+restjumpx;
                nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                nj(2)=njumpx-nj(1);
                %
                restjgrw=njumpx/2+restjgrw; % ~unbiased GRW
                njgrw(1)=floor(restjgrw); restjgrw=restjgrw-njgrw(1);
                njgrw(2)=njumpx-njgrw(1);
                %
                if x==1
                    nn(y,2)=nn(y,2)+nj(2);
                elseif x==I
                    nn(y,I-1)=nn(y,I-1)+nj(1);
                else
                    for i=1:2
                        xd=x+(2*i-3);
                        nn(y,xd)=nn(y,xd)+nj(i);

                        for m=1:length(mx)
                            if abs(x-imx(m))<=ia && abs(y-jmy)<ia && abs(xd-imx(m))<ia
                                dmxx(m)=dmxx(m)+nj(i)*xd*dx/(2*a)^2;
                                dmyy(m)=dmyy(m)+nj(i)*y*dy/(2*a)^2;
                                dmVx(m)=dmVx(m)+nj(i)*(xd-x)*dx/dt/(2*a)^2;
                                if normVx==0
                                    dmyVx(m)=dmyVx(m)+nj(i)*y*dy*(xd-x)*dx/dt/(2*a)^2;
                                    dmxVx(m)=dmxVx(m)+nj(i)*(x+xd)/2*dx*(xd-x)*dx/dt/(2*a)^2;
                                else
                                    dmyVx(m)=dmyVx(m)+njgrw(i)*y*dy*(xd-x)*dx/dt/(2*a)^2;
                                    dmxVx(m)=dmxVx(m)+njgrw(i)*(x+xd)/2*dx*(xd-x)*dx/dt/(2*a)^2;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%% Boundary conditions - concentration
%%%% BCY Left/Right
nn(:,1)=nn(:,2); % no flux left
nn(:,I)=nn(:,I-1); % no flux right
%%%% BCX Bottom/Top
nn(1,:)=nn(2,:); % no flux on bottom boundary
nn(J,:)=nn(J-1,:);  % no flux on top boundary

n=nn; c=n;

