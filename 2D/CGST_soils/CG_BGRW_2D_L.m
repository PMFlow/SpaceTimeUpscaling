% BGRW function - transport step in coupled flow & trnasport 2D

function [c]=CG_BGRW_2D_L(c0,cBC,dthtc,I,J,dx,dy,dtt,Vx,Vy,D1,D2,Lc,dJ,i1,i2)

%% BGRW solution
v=Vy*dtt/dx/Lc; u=Vx*dtt/dy/Lc;
ru=2*D1*dtt/Lc/dx^2*ones(J,I);
rv=2*D2*dtt/Lc/dy^2*ones(J,I);
n=c0; nBC=cBC;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0; restf=0;
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
                if(njumpy)>0
                    restjumpy=njumpy*0.5*(1-v(y,x)/ry)+restjumpy;
                    nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                    nj(2)=njumpy-nj(1);
                    if y==1
                        nn(2,x)=nn(2,x)+nj(2); 
                    elseif y==J
                        nn(J-1,x)=nn(J-1,x)+nj(1); 
                    else
                        for i=1:2
                            yd=y+(2*i-3);
                            nn(yd,x)=nn(yd,x)+nj(i);
                        end
                    end
                end                
                if(njumpx)>0
                    restjumpx=njumpx*0.5*(1-u(y,x)/rx)+restjumpx;
                    nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                    nj(2)=njumpx-nj(1);
                    if x==1
                        nn(y,2)=nn(y,2)+nj(2); 
                    elseif x==I
                        nn(y,I-1)=nn(y,I-1)+nj(1); 
                    else
                        for i=1:2
                            xd=x+(2*i-3);
                            nn(y,xd)=nn(y,xd)+nj(i);
                        end
                    end
                end
            end
        end
    end
    %% Boundary conditions - concentration
        %%%% BCY Left/Right
        nn(:,1)=nn(:,2); % no flux left % nBC(:,1); %
        nn(:,I)=nn(:,I-1); % no flux right % nBC(:,I); %
        %%%% BCX Bottom/Top
        nn(1,2:I-1)=nn(2,2:I-1); % nBC(1,2:I-1); % no flux on bottom boundary
        nn(J-dJ:J,i1:i2)=nBC(J-dJ:J,i1:i2); % set to initial condition     
    %% Source term concentration
    restf=dthtc(2:J-1,2:I-1)+restf; nf=floor(restf); restf=restf-nf;
    nn(2:J-1,2:I-1)=nn(2:J-1,2:I-1)+nf;
    n=nn; c=n;

