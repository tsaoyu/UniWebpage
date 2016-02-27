

Time_Test=xlsread('满载高速实船试验数据','A2:A58'); 
Veo_Ship_Test=xlsread('满载高速实船试验数据','B2:B58'); 
Tension1_Test=xlsread('满载高速实船试验数据','H2:H58'); 
Tension2_Test=xlsread('满载高速实船试验数据','I2:I58'); 
Tension_Test=Tension1_Test+Tension2_Test; 

%% 画图
figure; plot(Time_Test,Veo_Ship_Test,'k'); 
hold on; plot(t, Bv(1,:,N),'r')

figure; plot(Time_Test,Tension_Test,'k');  
hold on;  plot(t,sqrt(F_Hin_Kane(1,1:length(t),2).^2+F_Hin_Kane(3,1:length(t),2).^2),'r')