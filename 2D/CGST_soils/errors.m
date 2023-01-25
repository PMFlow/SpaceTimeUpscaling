
clear all

err2Dinf=zeros(2,3); err2D2=zeros(2,3); err2Dm=zeros(2,3);

load 2D_T132_centr.mat
disp('2 dim');
for i=1:3
    err2Dinf(1,i)=norm(mc1v(i,:)-mc1(i,:),inf)/norm(mc1(i,:),inf);
    err2Dinf(2,i)=norm(mc2v(i,:)-mc2(i,:),inf)/norm(mc2(i,:),inf);
    err2D2(1,i)=norm(mc1v(i,:)-mc1(i,:))/norm(mc1(i,:));
    err2D2(2,i)=norm(mc2v(i,:)-mc2(i,:))/norm(mc2(i,:));
    [~,I]=max(abs(mc1v(i,:)-mc1(i,:)));
    err2Dm(1,i)=max(abs(mc1v(i,:)-mc1(i,:)))/mc1(i,I);
    [~,I]=max(abs(mc2v(i,:)-mc2(i,:)));
    err2Dm(2,i)=max(abs(mc2v(i,:)-mc2(i,:)))/mc2(i,I);
end
disp('inf-norm')
disp(err2Dinf);
disp('2-norm')
disp(err2D2);
disp('max-error')
disp(err2Dm);

load 2D_T132_dec.mat
disp('2 dim-decentered');
for i=1:3
    err2Dinf(1,i)=norm(mc1v(i,:)-mc1(i,:),inf)/norm(mc1(i,:),inf);
    err2Dinf(2,i)=norm(mc2v(i,:)-mc2(i,:),inf)/norm(mc2(i,:),inf);
    err2D2(1,i)=norm(mc1v(i,:)-mc1(i,:))/norm(mc1(i,:));
    err2D2(2,i)=norm(mc2v(i,:)-mc2(i,:))/norm(mc2(i,:));
    [~,I]=max(abs(mc1v(i,:)-mc1(i,:)));
    err2Dm(1,i)=max(abs(mc1v(i,:)-mc1(i,:)))/mc1(i,I);
    [~,I]=max(abs(mc2v(i,:)-mc2(i,:)));
    err2Dm(2,i)=max(abs(mc2v(i,:)-mc2(i,:)))/mc2(i,I);
end
disp('inf-norm')
disp(err2Dinf);
disp('2-norm')
disp(err2D2);
disp('max-error')
disp(err2Dm);
