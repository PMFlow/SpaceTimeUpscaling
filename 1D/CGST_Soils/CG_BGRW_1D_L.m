function [c]=CG_BGRW_1D_L(c0,cBC,dthtc,I,dI,dx,dt,q,D,Lc)

%% BGRW solution
restr=0; restjump=0; restf=0; 
u=q*dt/dx;
ru=2*D*dt/Lc/dx^2*ones(1,I);
n=c0; nBC=cBC;
nn=zeros(1,I);

for x=1:I
    if n(x) > 0
        r=ru(x);
        restr=n(x)*(1-r)+restr; nsta=floor(restr);
        restr=restr-nsta; njump=n(x)-nsta;
        nn(x)=nn(x)+nsta;
        if(njump)>0
            restjump=njump*0.5*(1-u(x)/r)+restjump;
            nj(1)=floor(restjump); restjump=njump-nj(1);
            nj(2)=floor(restjump); restjump=restjump-nj(2);
            if x==1
                nn(2)=nn(2)+nj(2);
            elseif x==I
                nn(I-1)=nn(I-1)+nj(1);
            else
                for i=1:2
                    xd=x+(2*i-3);
                    nn(xd)=nn(xd)+nj(i);
                end
            end
        end
    end
end
%% boundary conditions
nn(1)=nn(2); % no flux % 0; % nBC(1); % 
% nn(I)=nBC(I); % set to IC 
nn(I-dI:I)=nBC(I-dI:I); % DEEP BC
%% Source term concentration
restf=dthtc(2:I-1)+restf; nf=floor(restf); restf=restf-nf;
nn(2:I-1)=nn(2:I-1)+nf;
c=nn;