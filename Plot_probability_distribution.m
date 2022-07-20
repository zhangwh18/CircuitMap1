clc
clear
data=load('/home/pcs/Documents/CircuitMap/TmeanA_0.500_2242716');
f_s5=data(:,1);
len=floor(1.1*sqrt(length(f_s5)));%calculate the number of sections
x=linspace(min(f_s5),max(f_s5),len);  %  divide the difference value into len sections
prob=histc(f_s5,x)/length(f_s5)/(x(2)-x(1));% find the probability into different sections 
% x=linspace(0,max(f_s5),15);%  divide the difference value into len sections
stairs(x,prob,'r-');

hold on
f_s5=data(:,2);
len=floor(1.1*sqrt(length(f_s5)));%calculate the number of sections
x=linspace(min(f_s5),max(f_s5),len);  %  divide the difference value into len sections
prob=histc(f_s5,x)/length(f_s5)/(x(2)-x(1));% find the probability into different sections 
% x=linspace(0,max(f_s5),15);%  divide the difference value into len sections
stairs(x,prob,'b-');
% data-data1;
% find(abs(ans)>0)
hold on
% data1=data
% endtime=700;
% time=data(1:endtime,1);
% Variance=data(1:endtime,2);
% loglog(time,Variance);
% sum((data(:,1)*1784391).^2+(data(:,2)*1784391).^2)