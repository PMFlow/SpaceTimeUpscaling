function [c1,c2,dmc1,dmc2,dmc1v,dmc2v]=CG_reaction_Monod2(c01,c02,dt,tht,Lc,I,J,my,imx,jmy,a,ia)

N=10^24;
M1=0.1; M2=0.1; a1=5; a2=0.5;  % ~ [Bause & Knabner, 2004]
c01=c01/N; c02=c02/N;
m=(c01./(M1+c01)).*(c02./(M2+c02));
c1=c01-dt*a1*m.*tht/Lc;
c2=c02-dt*a2*m.*tht/Lc;
c1=c1*N; c2=c2*N;

dmc1=zeros(1,length(my)); dmc2=zeros(1,length(my));
% dmc1p=zeros(1,length(my)); dmc2p=zeros(1,length(my));
dmc1v=zeros(1,length(my)); dmc2v=zeros(1,length(my));

for m=1:length(my)
    mi=1;
    for j=1:J
        for i=1:I
            if abs(i-imx)<ia && abs(j-jmy(m))<ia
                dmc1(m)=dmc1(m)+c1(j,i)/(2*a)^2; dmc1v(m)=dmc1(m);
                dmc2(m)=dmc2(m)+c2(j,i)/(2*a)^2; dmc2v(m)=dmc2(m);
                mi=mi+1;
            end
            %         dmc1p(m)=dmc1(m)*(2*a)/mi; dmc2p(m)=dmc2(m)*(2*a)/mi;
        end
    end
    %     fprintf('mi =  %d \n',mi) ;
end
