%% Computation of maximum relative differences between volume and CGST averages

clear all

err1Dinf=zeros(2,3); err1D2=zeros(2,3); err1Dm=zeros(2,3);

load CG_aquifer.mat
disp('1 dim-aquifer');
for i=1:3
    err1Dinf(1,i)=norm(mc1v(i,:)-mc1(i,:),inf)/norm(mc1(i,:),inf);
    err1Dinf(2,i)=norm(mc2v(i,:)-mc2(i,:),inf)/norm(mc2(i,:),inf);
    err1D2(1,i)=norm(mc1v(i,:)-mc1(i,:))/norm(mc1(i,:));
    err1D2(2,i)=norm(mc2v(i,:)-mc2(i,:))/norm(mc2(i,:));
    [~,I]=max(abs(mc1v(i,:)-mc1(i,:)));
    err1Dm(1,i)=max(abs(mc1v(i,:)-mc1(i,:)))/mc1(i,I);
    [~,I]=max(abs(mc2v(i,:)-mc2(i,:)));
    err1Dm(2,i)=max(abs(mc2v(i,:)-mc2(i,:)))/mc2(i,I);
end
disp('inf-norm')
disp(err1Dinf);
disp('2-norm')
disp(err1D2);
disp('max-error')
disp(err1Dm);

load CG_aquifer_ens.mat
disp('1 dim-aquifer_ensemble');
for i=1:3
    err1Dinf(1,i)=norm(mc1v(i,:)/R-mc1(i,:)/R,inf)/norm(mc1(i,:)/R,inf);
    err1Dinf(2,i)=norm(mc2v(i,:)/R-mc2(i,:)/R,inf)/norm(mc2(i,:)/R,inf);
    err1D2(1,i)=norm(mc1v(i,:)/R-mc1(i,:)/R)/norm(mc1(i,:)/R);
    err1D2(2,i)=norm(mc2v(i,:)/R-mc2(i,:)/R)/norm(mc2(i,:)/R);
    [~,I]=max(abs(mc1v(i,:)/R-mc1(i,:)/R));
    err1Dm(1,i)=max(abs(mc1v(i,:)/R-mc1(i,:)/R))/(mc1(i,I)/R);
    [~,I]=max(abs(mc2v(i,:)/R-mc2(i,:)/R));
    err1Dm(2,i)=max(abs(mc2v(i,:)/R-mc2(i,:)/R))/(mc2(i,I)/R);
end
disp('inf-norm')
disp(err1Dinf);
disp('2-norm')
disp(err1D2);
disp('max-error')
disp(err1Dm);
