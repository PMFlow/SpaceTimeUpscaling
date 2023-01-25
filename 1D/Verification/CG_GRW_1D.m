function [c,dmxx,dmvv,dmxv]=CG_GRW_1D(c0,I,dx,dt,q,D,mx,imx,a,ia,d,stepU)

%% (unbiased) GRW solution
u=floor(q*stepU+0.5);
r=2*D*dt/d^2/dx^2*ones(1,I);
n=c0;
nn=zeros(1,I);
restr=0; restjump=0;
dmxx=zeros(1,length(mx)); dmvv=zeros(1,length(mx)); dmxv=zeros(1,length(mx));

for x=1:I
    if n(x) > 0
        xa=x+u(x);
        restr=n(x)*(1-r(x))+restr; nsta=floor(restr);
        restr=restr-nsta; njump=n(x)-nsta;
        if xa<1 
            nn(1)=nn(1)+n(x);
        elseif xa>I
            nn(I)=nn(I)+n(x);
        else
            nn(xa)=nn(xa)+nsta;
            for m=1:length(mx)
                if abs(x-imx(m))<=ia && abs(xa-imx(m))<ia    
                    dmxx(m)=dmxx(m)+n(x)*xa*dx/(2*a);
                    dmvv(m)=dmvv(m)+n(x)*(xa-x)*dx/dt/(2*a); % translation with constant velocity; unbiased jumps do not contribute to <\xi> !                 
                end
            end
            if(njump)>0
                restjump=njump/2+restjump;
                nj(1)=floor(restjump); restjump=restjump-nj(1);
                nj(2)=njump-nj(1);
                if xa==1
                    nn(2)=nn(2)+nj(2);
                elseif xa==I
                    nn(I-1)=nn(I-1)+nj(1);
                else
                    for i=1:2
                        if nj(i)>0
                            xj=xa+(2*i-3)*d;
                            if xa~=x && xj<=1
                                nn(xa)=nn(xa)+nj(2);
                            elseif xa~=x && xj>I
                                nn(xa)=nn(xa)+nj(1);
                            else
                                nn(xj)=nn(xj)+nj(i);
                                for m=1:length(mx)
                                    if abs(xa-imx(m))<=ia && abs(xj-imx(m))<ia
%                                         dmxx(m)=dmxx(m)+nj(i)*xj*dx/(2*a);
                                        dmxv(m)=dmxv(m)+nj(i)*(xa+xj)/2*dx*(xj-xa)*dx/dt/(2*a);
                                    end
                                end
                            end %
                        end
                    end
                end
            end
        end
    end
end
%% boundary conditions
nn(1)=nn(2); % no flux
nn(I)=nn(I-1);
c=nn;
