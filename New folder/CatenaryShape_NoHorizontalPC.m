
% CatenaryShape_NoHorizontalPC为求解无卧链悬链线平面形状的程序
% 已知缆索长L_Pre、缆索上端至海底距离h_Cha、缆索水平投影x_Pre、锚链水中单位重量w_PUW，求解其他物理量；
% 注意坐标系的问题：本程序的坐标系与主程序DynRes_AncCha_Anc1中的惯性系不同
% 本程序中坐标系原点位于左下端缆索处，求解缆索力及形状参数等；
% 坐标的变换在主程序中进行☆☆☆

global  L_Pre  h_Pre  x_Pre  x_PreS  w_PUW  N_Pre  X_HoP  Z_HoP  Tv_Pre  theta_Pre

%% 基于fsolve函数求解联立非线性方程组(1)(2),solve不能找到显式解法
% (1) a*cosh(x/a+log(tan(theta0)+sec(theta0)))-a*sec(theta0)-z
% (2) a*sinh(x/a+log(tan(theta0)+sec(theta0)))-a*tan(theta0)-L_Pre
% (3) Th*cosh(x/a+log(tan(theta0)+sec(theta0)))-T
% 对xS重新赋值
x_PreS=x_Pre; 
% ★注意求解初值的选择十分重要，不合适的初值可能会求解错误
sat0=[10, 0];
[sat,fval,exitflag]=fsolve(@ NoHorizontalPC,sat0);
a=sat(1);  theta0=sat(2);   
Th_Pre=a*w_PUW;
T=Th_Pre*cosh(x_PreS/a+log(tan(theta0)+sec(theta0)));
Tv_Pre=sqrt(T^2-Th_Pre^2);

%% 计算缆索各铰接点的x,z坐标 
%  X_HoP, Z_HoP分别为缆索各铰链点的x,z坐标矩阵(列向量)
X_HoP=zeros(N_Pre+1,1);  Z_HoP=zeros(N_Pre+1,1);
for k=1:N_Pre+1
    % 通过s=0～L_Pre，求解锚链各个铰接点的坐标
    % ★注意：锚链各段的编号从左下端开始1,2,...★
    s=(k-1)*L_Pre/N_Pre;
    x_int=fsolve(@(x)a*sinh(x/a+log(tan(theta0)+sec(theta0)))-a*tan(theta0)-s,x_PreS);
    z_int=a*cosh(x_int/a+log(tan(theta0)+sec(theta0)))-a*sec(theta0);
    X_HoP(k)=x_int;
    Z_HoP(k)=z_int;
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
% title('No Horizontal Prevention Cable');  
