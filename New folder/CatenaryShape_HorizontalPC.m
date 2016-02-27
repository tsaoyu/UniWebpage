
% CatenaryShape_HorizontalPC为求解有卧链悬链线平面形状的程序
% 已知缆索长度L_Pre、缆索上端至海底距离h_Pre、缆索水平投影x_Pre、缆索水中单位重量w_PUW，求解其他物理量；
% 注意坐标系的问题：本程序的坐标系与主程序DynRes_AncCha_Anc1中的惯性系不同
% 本程序中坐标系原点位于左下端缆索处，求解缆索力及形状参数等；
% 坐标的变换在主程序中进行☆☆☆

global  L_Pre  h_Pre  x_Pre  x_PreS  w_PUW  N_Pre  X_HoP  Z_HoP  Tv_Pre  theta_Pre

%% 基于fsolve函数求解其他参数（形状参数，拉力等），
% 先按坐标系原点位于左下端锚链求解形状参数，拉力等 
% b0为b的初值
b0=10;       
b=fsolve(@(b)(L_Pre-x_Pre)*(cosh(b)-1)-h_Pre*(sinh(b)-b),b0);
% 计算其他物理量 
% a=Th_Pre_Pre/w_PUW, b=x_PreS/a均为便于求解而设置的参数; 
Th_Pre=w_PUW*(L_Pre-x_Pre)/(sinh(b)-b); 
Tv_Pre=Th_Pre*sinh(b); 
a=Th_Pre/w_PUW;  
Sc=sinh(b)*a;
x_PreS=a*b;

%% 计算锚链各铰接点的x,z坐标 
%  X_HoP, Z_HoP分别为锚链各铰链点的x,z坐标矩阵(列向量)
X_HoP=zeros(N_Pre+1,1);  Z_HoP=zeros(N_Pre+1,1);
for k=1:N_Pre+1
    % 通过s=0～L_Pre，求解锚链各个铰接点的坐标
    % ★注意：锚链各段的编号从左下端开始1,2,...★
    s=(k-1)*L_Pre/N_Pre;
    if s<=x_Pre-x_PreS
        X_HoP(k)=s;           Z_HoP(k)=0;
    else
        Sc_s=s-(L_Pre-Sc);            xS_s=a*asinh(Sc_s/a);
        X_HoP(k)=x_Pre-(x_PreS-xS_s);        Z_HoP(k)=a*(cosh(xS_s/a)-1);
    end
end
% 计算各链环初始角度(列向量)
theta_Pre=zeros(N_Pre,1);
for k=1:N_Pre
    theta_Pre(k)=atan((Z_HoP(k+1)-Z_HoP(k))/(X_HoP(k+1)-X_HoP(k)));
    if abs(theta_Pre(k)-pi/2)<1.0e-6||abs(theta_Pre(k)+pi/2)<1.0e-6
        pause;   disp('错误！缆索出现竖直段！')
    end
    %     % 转换为与x轴正向的钝角(负值)
end
% set(0,'DefaultFigureColor','w');
% set(0,'DefaultAxesFontname','Times New Roman');
% figure;  plot(X_HoP,Z_HoP,'-o','LineWidth',1); 
% title('Exist Horizontal Prevention Cable');