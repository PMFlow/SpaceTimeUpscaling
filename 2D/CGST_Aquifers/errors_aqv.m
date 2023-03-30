%% Computation of relative differences between volume and CGST averages
function errors_aqv(mc1,mc1v,mc2,mc2v)

err1Dinf=zeros(2,3); err1D2=zeros(2,3); err1Dm=zeros(2,3);
disp('2 dim-aquifer');
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
disp('l2-norm')
disp(err1D2);
disp('max-error')
disp(err1Dm);
