clc
clear
data=load('/home/pcs/Documents/CircuitMap/Variance_g_1_N_5000_M_1_T_7');
endtime=700;
time=data(1:endtime,1);
Variance=data(1:endtime,2);
loglog(time,Variance);
hold on
data=load('/home/pcs/Documents/CircuitMap/Variance_g_1_N_5000_M_1_T_5_s_05');
time=data(1:endtime,1);
Variance=data(1:endtime,2);
loglog(time,Variance);
% data=load('/home/pcs/Documents/CircuitMap/HFt0_try_ini');
% time=data(1:endtime,1);
% Variance=data(1:endtime,2);
% loglog(time,Variance);
% hold on
time1=time(9:280);
loglog(time1,1*time1.^(-0.5),'--')
time2=time(140:588);
% loglog(time2,0.005*time2.^(-1.8),'--')

% legend('0.5','0.05','0','T^{-0.5}','T^{-1.8}')
legend('0.5','T^{-0.5}')
xlabel('T')
ylabel('\sigma^2(LO)')