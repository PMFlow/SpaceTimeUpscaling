function [c,dmxx,dmyy,dmVx,dmVy,dmxVx,dmxVy,dmyVx,dmyVy]=CG_GRW_2D(c0,I,J,dx,dy,dt,Vx,Vy,D1,D2,mx,a,ia,imx,jmy,d,stepU)

%% GRW solution
vr=Vy*dt/dx; ur=Vx*dt/dy;
u=floor(ur*stepU+0.5);
v=floor(vr*stepU+0.5);
rx=2*D1*dt/(d^2*dx^2); ry=2*D2*dt/(d^2*dy^2); r=rx+ry;
n=c0;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0;
dmxx=zeros(1,length(mx)); dmyy=zeros(1,length(mx));
dmVx=zeros(1,length(mx)); dmVy=zeros(1,length(mx));
dmxVx=zeros(1,length(mx)); dmxVy=zeros(1,length(mx));
dmyVx=zeros(1,length(mx)); dmyVy=zeros(1,length(mx));
if max(max(r))>1
    disp(r);
    return
end
for y=1:J
    for x=1:I

        if n(y,x) > 0
            xa=x+u(y,x); ya=y+v(y,x);
            restr=n(y,x)*(1-r)+restr; nsta=floor(restr);
            restr=restr-nsta; njump=n(y,x)-nsta;
            if ya<1 || ya>J
                ya=y;
            end
            if xa<1 || xa>I
                xa=x;
            end
            if ya<1
                ya=1;
            end
            if ya>J
                ya=J;
            end
            if xa<1
                xa=1;
            end
            if xa>I
                xa=I;
            end
            nn(ya,xa)=nn(ya,xa)+nsta;
            restjump=njump*rx/r+restjump;
            njumpx=floor(restjump); restjump=restjump-njumpx;
            njumpy=njump-njumpx;
            for m=1:length(mx)
                if abs(x-imx(m))<=ia && abs(xa-imx(m))<ia && abs(y-jmy)<=ia && abs (ya-jmy)<ia
                    dmxx(m)=dmxx(m)+nsta*xa*dx/(2*a)^2;
                    dmyy(m)=dmyy(m)+nsta*ya*dy/(2*a)^2;
                    dmVy(m)=dmVy(m)+n(y,x)*(ya-y)*dy/dt/(2*a)^2; % translation with constant velocity; unbiased jumps do not contribute to <\xi> !
                    dmVx(m)=dmVx(m)+n(y,x)*(xa-x)*dx/dt/(2*a)^2;
                end

            end
            if(njumpy)>0
                restjumpy=njumpy/2+restjumpy;
                nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                nj(2)=njumpy-nj(1);
                if ya==1
                    nn(2,xa)=nn(2,xa)+nj(2);
                elseif ya==J
                    nn(J-1,xa)=nn(J-1,xa)+nj(1);
                else
                    for i=1:2
                        if nj(i)>0
                            yd=ya+(2*i-3)*d;
                            if ya~=y && yd<=1
                                nn(ya,xa)=nn(ya,xa)+nj(2);
                            elseif ya~=y && yd>J
                                nn(ya,xa)=nn(ya,xa)+nj(1);
                            else
                                nn(yd,xa)=nn(yd,xa)+nj(i);
                                for m=1:length(mx)
                                    if abs(xa-imx(m))<ia && abs(ya-jmy)<=ia && abs(yd-jmy)<ia
                                        dmyy(m)=dmyy(m)+nj(i)*yd*dy/(2*a)^2;
                                        dmxx(m)=dmxx(m)+nj(i)*x*dx/(2*a)^2;
                                        dmyVy(m)=dmyVy(m)+nj(i)*(ya+yd)/2*dy*(yd-ya)*dy/dt/(2*a)^2;
                                        dmxVy(m)=dmxVy(m)+nj(i)*xa*dx*(yd-ya)*dy/dt/(2*a)^2;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if(njumpx)>0
                restjumpx=njumpx/2+restjumpx;
                nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                nj(2)=njumpx-nj(1);
                if xa==1
                    nn(ya,2)=nn(ya,2)+nj(2);
                elseif xa==I
                    nn(ya,I-1)=nn(ya,I-1)+nj(1);
                else
                    for i=1:2
                        if nj(i)>0
                            xd=xa+(2*i-3)*d;
                            if xa~=x && xd<=1
                                nn(ya,xa)=nn(ya,xa)+nj(2);
                            elseif xa~=x && xd>I
                                nn(ya,xa)=nn(ya,xa)+nj(1);
                            else
                                nn(ya,xd)=nn(ya,xd)+nj(i);
                                for m=1:length(mx)
                                    if abs(xa-imx(m))<=ia && abs(ya-jmy)<ia && abs(xd-imx(m))<ia
                                        dmxx(m)=dmxx(m)+nj(i)*xd*dx/(2*a)^2;
                                        dmyy(m)=dmyy(m)+nj(i)*y*dy/(2*a)^2;
                                        dmyVx(m)=dmyVx(m)+nj(i)*ya*dy*(xd-xa)*dy/dt/(2*a)^2;
                                        dmxVx(m)=dmxVx(m)+nj(i)*(xa+xd)/2*dx*(xd-xa)*dy/dt/(2*a)^2;
                                    end
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

n=nn;
c=n;
