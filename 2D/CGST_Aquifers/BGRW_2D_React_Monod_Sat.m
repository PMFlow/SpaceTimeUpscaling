%% BGRW function (random) - transport step in coupled flow & trnasport 2D
function [c]=BGRW_2D_React_Monod_Sat(c0,cBC,tht,I,J,dx,dy,i0,j0,dtt,Vx,Vy,D1,D2,dw)

%% BGRW solution
v=Vy*dtt/dx./tht; u=Vx*dtt/dy./tht;
ru=2*D1*dtt/dx^2*ones(J,I)./tht; 
rv=2*D2*dtt/dy^2*ones(J,I)./tht;
n=c0; nBC=cBC;
nn=zeros(J,I);
restr=0; restjump=0; restjumpx=0; restjumpy=0; 
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
%%%% BCXBottom/Top
nn(:,1:j0-dw-1)=nBC(:,1:j0-dw-1); % set to IC on bottom boundary %% J=>x; I=>y !!!
nn(:,I-j0+dw+1:I)=nBC(:,I-j0+dw+1:I); % set to IC on top boundary
%%%% BCYLeft/Right
nn(1:j0-dw-1,2:I-1)=nBC(1:j0-dw-1,2:I-1);  %  set to IC on left boundary
nn(j0-dw:j0+dw,i0-dw:i0+dw)=nBC(j0-dw:j0+dw,i0-dw:i0+dw);  % injection
nn(J,2:I-1)=nn(J-1,2:I-1); % no flux on rigth boundary %  
n=nn; c=n;
