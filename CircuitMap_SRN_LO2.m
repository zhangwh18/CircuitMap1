clc
clear
N=10^4/2;

% theta=0;
g=1;
T=10^9;
M=2000;
time=linspace(0,10^9,10^9+1);
% thetas=[0.5,0.25,0.1,0.075,0.05,0.025,0.01,0.0075,0.005,0.0025,0.001];
thetas=[0.5,0.05,0];

% time=[10^0:10^0:10^1,10^1:10^1:10^2,10^2:10^2:10^3,10^3:10^3:10^4];
for theta=thetas
            x0A=rand(M,N)+1i.*rand(M,N);
            x0B=rand(M,N)+1i.*rand(M,N);

            for i=1:M
                
            x0A(M,:)=x0A(M,:)/sum(abs(x0A(M,:)));
                        x0B(M,:)=x0B(M,:)/sum(abs(x0B(M,:)));
            end
phiA=zeros(M,N)+1i.*zeros(M,N);
phiB=zeros(M,N)+1i.*zeros(M,N);
MeanvalueA=abs(x0A);
%     x0A=zeros(1,N)+1i.*zeros(1,N);
%     x0B=zeros(1,N)+1i.*zeros(1,N);

%     x0A(1)=1+0*1i;

% x0A(1)=rand(1,N)+1i*;
% Meanvalue=zeros(1,5000);
Variance=zeros(1,length(time));
timeind=0;
for t=time

    timeind=timeind+1;
%     for M=1:100
%     Meanvalue(M)=mean(abs(x0A));

%     Variance(timeind)=mean((abs(x0A)-Meanvalue(timeind)).^2);
    
    for n=1:N
        if n==1
            phiA(:,n)=cos(theta)^2*x0A(:,n)-cos(theta)*sin(theta)*x0B(:,N)+sin(theta)^2*x0A(:,n+1)+cos(theta)*sin(theta)*x0B(:,n);
            phiB(:,n)=cos(theta)^2*x0B(:,n)-cos(theta)*sin(theta)*x0A(:,n)+sin(theta)^2*x0B(:,N)+cos(theta)*sin(theta)*x0A(:,n+1);
        elseif n==N
            phiA(:,n)=cos(theta)^2*x0A(:,n)-cos(theta)*sin(theta)*x0B(:,n-1)+sin(theta)^2*x0A(:,1)+cos(theta)*sin(theta)*x0B(:,n);
            phiB(:,n)=cos(theta)^2*x0B(:,n)-cos(theta)*sin(theta)*x0A(:,n)+sin(theta)^2*x0B(:,n-1)+cos(theta)*sin(theta)*x0A(:,1);
            
        else
            phiA(:,n)=cos(theta)^2*x0A(:,n)-cos(theta)*sin(theta)*x0B(:,n-1)+sin(theta)^2*x0A(:,n+1)+cos(theta)*sin(theta)*x0B(:,n);
            phiB(:,n)=cos(theta)^2*x0B(:,n)-cos(theta)*sin(theta)*x0A(:,n)+sin(theta)^2*x0B(:,n-1)+cos(theta)*sin(theta)*x0A(:,n+1);
        end
    end
    x0A=exp(1i.*g.*abs(phiA).^2).*phiA;
    x0B=exp(1i.*g.*abs(phiB).^2).*phiB;
%     end
    
    MeanvalueA=(abs(MeanvalueA)*timeind+abs(x0A))/(timeind+1);
    
    Variance(timeind)=mean((abs(MeanvalueA(:,2))-mean(abs(MeanvalueA(:,2)))).^2);
%     Variance(timeind)=mean((Meanvalue-mean(Meanvalue))).^2;
end
figure(2)
loglog(time,Variance);
save(['/home/pcs/Documents/CircuitMap/data/variance_',num2str(theta),num2str(M),num2str(N),'.mat'],'Variance','-v7.3');
hold on
% legend(num2str(theta,4))
end
% legend([num2str(thetas,4)])
legend('0.5','0.001','0')
% legend('0.5','0.25','0.1','0.075','0.05','0.025','0.01','0.0075','0.005','0.0025','0.001')