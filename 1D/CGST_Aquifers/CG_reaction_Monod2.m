function [c1,c2,dmc1,dmc2,dmc1p,dmc2p,dmc1v,dmc2v]=CG_reaction_Monod2(c01,c02,dt,I,mx,imx,a,ia,NI1,NI2)

K1=0.1; K2=0.1; Kr1=5; Kr2=0.5;  % ~ [Bause & Knabner, 2004]

c01=c01/NI1; c02=c02/NI2;
c1(1)=c01(1); c1(I)=c01(I);
c2(1)=c02(1); c2(I)=c02(I);

c_bio=1;
m=c_bio*(c01./(K1+c01)).*(c02./(K2+c02));
c1(2:I-1)=c01(2:I-1)-dt*Kr1*m(2:I-1);
c2(2:I-1)=c02(2:I-1)-dt*Kr2*m(2:I-1);

c1=c1*NI1; c2=c2*NI2;

dmc1=zeros(1,length(mx)); dmc2=zeros(1,length(mx));
dmc1p=zeros(1,length(mx)); dmc2p=zeros(1,length(mx));
dmc1v=zeros(1,length(mx)); dmc2v=zeros(1,length(mx)); 

for m=1:length(mx)
    mi=1;
    for i=1:I
        if abs(i-imx(m))<ia
            dmc1(m)=dmc1(m)+c1(i)/(2*a); dmc1v(m)=dmc1(m);
            dmc2(m)=dmc2(m)+c2(i)/(2*a); dmc2v(m)=dmc2(m);
            mi=mi+1;
        end
        dmc1p(m)=dmc1(m)*(2*a)/mi; dmc2p(m)=dmc2(m)*(2*a)/mi;
    end
%     fprintf('mi =  %d \n',mi) ;
end
