clc
clear
N=10^5/2;

% theta=0;
theta=0.33*pi;
T=10^9;
time=logspace(1,9,1000);
gs=[0.5,0.1,0.005,0.0025];
% time=[10^0:10^0:10^1,10^1:10^1:10^2,10^2:10^2:10^3,10^3:10^3:10^4];
for g=gs
    x0A=zeros(1,N)+1i.*zeros(1,N);
    x0B=zeros(1,N)+1i.*zeros(1,N);
    phiA=zeros(1,N)+1i.*zeros(1,N);
    phiB=zeros(1,N)+1i.*zeros(1,N);
    x0A(1)=1+0*1i;
    Meanvalue=zeros(1,length(time));
    Variance=zeros(1,length(time));
    timeind=0;
    for t=time
        timeind=timeind+1;
        Meanvalue(timeind)=mean([abs(x0A),abs(x0B)]);
        Variance(timeind)=mean(([abs(x0A),abs(x0B)]-Meanvalue(timeind)).^2);
        
        for n=1:N
            if n==1
                phiA(n)=cos(theta)^2*x0A(n)-cos(theta)*sin(theta)*x0B(N)+sin(theta)^2*x0A(n+1)+cos(theta)*sin(theta)*x0B(n);
                phiB(n)=cos(theta)^2*x0B(n)-cos(theta)*sin(theta)*x0A(n)+sin(theta)^2*x0B(N)+cos(theta)*sin(theta)*x0A(n+1);
            elseif n==N
                phiA(n)=cos(theta)^2*x0A(n)-cos(theta)*sin(theta)*x0B(n-1)+sin(theta)^2*x0A(1)+cos(theta)*sin(theta)*x0B(n);
                phiB(n)=cos(theta)^2*x0B(n)-cos(theta)*sin(theta)*x0A(n)+sin(theta)^2*x0B(n-1)+cos(theta)*sin(theta)*x0A(1);
                
            else
                phiA(n)=cos(theta)^2*x0A(n)-cos(theta)*sin(theta)*x0B(n-1)+sin(theta)^2*x0A(n+1)+cos(theta)*sin(theta)*x0B(n);
                phiB(n)=cos(theta)^2*x0B(n)-cos(theta)*sin(theta)*x0A(n)+sin(theta)^2*x0B(n-1)+cos(theta)*sin(theta)*x0A(n+1);
            end
        end
        x0A=exp(1i.*g.*abs(phiA).^2).*phiA;
        x0B=exp(1i.*g.*abs(phiB).^2).*phiB;
    end
    loglog(time,abs(Variance)/max(abs(Variance)));
    hold on
    % legend(num2str(theta,4))
end
% legend([num2str(thetas,4)])
% legend('0.5','0.25','0.1','0.075','0.05','0.025','0.01','0.0075','0.005','0.0025','0.001')
legend('0.5','0.1','0.005','0.0025')