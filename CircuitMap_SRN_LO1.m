clc
clear
N=10^4/2;

% theta=0;
g=1;
T=10^9;
time=logspace(1,9,1000);
thetas=[0.5,0.25,0.1,0.075,0.05,0.025,0.01,0.0075,0.005,0.0025,0.001];
% time=[10^0:10^0:10^1,10^1:10^1:10^2,10^2:10^2:10^3,10^3:10^3:10^4];
for theta=thetas
    x0A=rand(1,N)+1i.*rand(1,N);
x0B=rand(1,N)+1i.*rand(1,N);
phiA=zeros(1,N)+1i.*zeros(1,N);
phiB=zeros(1,N)+1i.*zeros(1,N);
% x0A(1)=rand(1,N)+1i*;
Meanvalue=zeros(1,length(time));
Variance=zeros(1,length(time));
timeind=0;
for t=time
    timeind=timeind+1;
    Meanvalue(timeind)=mean(abs(x0A));
    Variance(timeind)=mean((abs(x0A)-Meanvalue(timeind)).^2);
    
    for n=1:N
        if n==1
            phiA(n)=x0A(n)-theta*(x0B(N)-x0B(n));
            phiB(n)=x0B(n)+theta*(x0A(n+1)-x0A(n)); 
        elseif n==N
            phiA(n)=x0A(n)-theta*(x0B(n-1)-x0B(n));
            phiB(n)=x0B(n)+theta*(x0A(1)-x0A(n)); 
            
        else
            phiA(n)=x0A(n)-theta*(x0B(n-1)-x0B(n));
            phiB(n)=x0B(n)+theta*(x0A(n+1)-x0A(n));        
        end

    end
    x0A=exp(1i.*g.*abs(phiA).^2).*phiA;
    x0B=exp(1i.*g.*abs(phiB).^2).*phiB;
end
loglog(time,Variance);
hold on
% legend(num2str(theta,4))
end
% legend([num2str(thetas,4)])
legend('0.5','0.25','0.1','0.075','0.05','0.025','0.01','0.0075','0.005','0.0025','0.001')