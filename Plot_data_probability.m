clc
clear
data=load('/home/pcs/Documents/CircuitMap/datac.tab');
x=data(:,1);
y=data(:,2);
phi=data(:,3);
f=x.^2.*exp(-x.*y.^2-y.^2+2.*y-4.*x);
% plot(x,f,'.')
f_s5=x;
len=floor(1.1*sqrt(length(f_s5)));%calculate the number of sections
x=linspace(min(f_s5),max(f_s5),len);  %  divide the difference value into len sections
prob=histc(f_s5,x)/length(f_s5)/(x(2)-x(1));% find the probability into different sections 
% x=linspace(0,max(f_s5),15);%  divide the difference value into len sections
stairs(x,prob,'-');
hold on 
plot(x,exp(-x));
endtime=588;