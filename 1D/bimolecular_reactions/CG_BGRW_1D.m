function [c,dmxx,dmvv,dmxv]=CG_BGRW_1D(c0,cBC,I,dx,dt,q,D,mx,a,imx,ia)

%% BGRW solution
restr=0; restjump=0; restjgrw=0; 
u=q*dt/dx;
ru=2*D*dt/dx^2*ones(1,I);
n=c0; nBC=cBC;
nn=zeros(1,I);
dmxx=zeros(1,length(mx)); dmvv=zeros(1,length(mx)); dmxv=zeros(1,length(mx));

for x=1:I
    if n(x) > 0
        r=ru(x);
        restr=n(x)*(1-r)+restr; nsta=floor(restr);
        restr=restr-nsta; njump=n(x)-nsta;
        nn(x)=nn(x)+nsta;
        for m=1:length(mx)
            if abs(x*dx-mx(m))<a
                dmxx(m)=dmxx(m)+nsta*x*dx/(2*a);
            end
        end
        if(njump)>0
            restjump=njump*0.5*(1-u(x)/r)+restjump;
            nj(1)=floor(restjump); restjump=njump-nj(1);
            nj(2)=floor(restjump); restjump=restjump-nj(2);
            %
            restjgrw=njump/2+restjgrw; % ~unbiased GRW
            njgrw(1)=floor(restjgrw); restjgrw=restjgrw-njgrw(1);
            njgrw(2)=njump-njgrw(1);
            %                        
            
            if x==1
                nn(2)=nn(2)+nj(2);
            elseif x==I
                nn(I-1)=nn(I-1)+nj(1);
            else
                for i=1:2
                    xd=x+(2*i-3);
                    nn(xd)=nn(xd)+nj(i);
                    for m=1:length(mx)
                        if abs(x-imx(m))<=ia && abs(xd-imx(m))<ia   
                            dmxx(m)=dmxx(m)+nj(i)*xd*dx/(2*a);
                            dmvv(m)=dmvv(m)+nj(i)*(xd-x)*dx/dt/(2*a);
                            dmxv(m)=dmxv(m)+(x+xd)/2*dx*njgrw(i)*(xd-x)*dx/dt/(2*a);
                        end
                    end
                end
            end
        end
    end
end
%% boundary conditions
nn(1)=nBC(1); % BC
nn(I)=nn(I-1); % no flux
c=nn;