function [c1,c2,dmc1,dmc2,dmc1p,dmc2p,dmc1v,dmc2v]=CG_reaction(c01,c02,dt,I,mx,imx,a,ia)

Kr1=0.1; Kr2=Kr1; 
N1=10^24; N2=N1;
c01=c01/N1; c02=c02/N2;

c1(1)=c01(1); c1(I)=c01(I);
c2(1)=c02(1); c2(I)=c02(I);

c1(2:I-1)=c01(2:I-1)-dt*Kr1*c01(2:I-1).*c02(2:I-1).^2;
c2(2:I-1)=c02(2:I-1)+dt*Kr2*c01(2:I-1).*c02(2:I-1).^2; 

c1=c1*N1; c2=c2*N2;

dmc1=zeros(1,length(mx));  dmc2=zeros(1,length(mx));
dmc1p=zeros(1,length(mx)); dmc2p=zeros(1,length(mx)); 
dmc1v=zeros(1,length(mx)); dmc2v=zeros(1,length(mx)); 

for m=1:length(mx)
    mi=0;
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
