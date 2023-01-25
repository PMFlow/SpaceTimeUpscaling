function [Vx,Vy] = velocity(I,J,dx,dy,p0,D)

%%       V_interior of \Omega        
        Vx(2:J-1,2:I-1)=-D.*((p0(2:J-1,3:I)-p0(2:J-1,1:I-2))/(2*dx)); 
        Vy(2:J-1,2:I-1)=-D.*((p0(3:J,2:I-1)-p0(1:J-2,2:I-1))/(2*dy)+1); 
%%       V_normal to boundary extended from interior of \Omega
        Vx(:,1)=Vx(:,2); Vx(:,I)=Vx(:,I-1); 
        Vx(1,:)=Vx(2,:); Vx(J,:)=Vx(J-1,:); 
        Vy(:,1)=Vy(:,2); Vy(:,I)=Vy(:,I-1);
        Vy(1,:)=Vy(2,:); Vy(J,:)=Vy(J-1,:); 
