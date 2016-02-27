
% 全套程序是基于Kane方法求解锚链的动力响应(2维，刚开始重力锚不动，然后走锚)
% 本主程序的主要功能是给出相关输入变量,求出积分函数的初值
% 然后基于龙格库塔法解一阶常微分方程组，最后进行动画演示

% 注意：在船舶与锚链连接之间加入一段竖直杆件以便定义船舶的初速度！☆☆☆
% 重力锚锚块被简化为长方体，长宽高分别为1.4m*5.6m*2.7m（单个12.5t为1.4m*1.4m*2.7m）
% 锚链为无挡锚链，单根链径为60mm,两个锚链的单位质量m_CU为2*71.756=143.5120kg/m, 湿重w_CUW=-m_CU*g-FB_CU  
% 锚链各节链环被简化为圆柱体杆件，等效半径R_CE=sqrt(m_CU/7.85e3/pi)=0.0763m  
% 拦阻索等效为钢丝绳，钢丝直径7mm，19根所以模型中的等效半径R_PE=sqrt(m_PU/7.85e3/pi)=0.0167； (2) 0.0035m*sqrt(38)=0.0216m
% 船舶被简化为长方体，长宽高分别为40m*10m*5m（500t时吃水1.25m,850t时吃水2.125m）

% 建立静止坐标系：固定点o0为原点，x0轴水平向右，z0轴竖直向上，y0轴竖直向内
% 静止坐标系原点位于船舶正浮状态的水线面形心处（可不与质心重合）
% 以海底面为重力势能的零势能面

clc;  clear ;  close all;

%% 定义全局变量
global Delta_t  Step_n  c_Seabed  g  rou_W  C_D  C_M  C_F  ks  upsilon  BvSFlag 
global N  Free  Low
global N_Cha  L_Pre  h_Pre  x_Pre  w_PUW  N_Pre  X_HoP  Z_HoP  Tv_Pre  theta  
global R_Cha  Len1  Len_Cha  Len_Pre  L_Aux  LenN  Len  Rad  Vol_Anc  M_Shi  M_Anc  Mass 
global Dis  S_Wet  IT_Local Q_Local  R_Local      
global F  M  F_Local  M_Local  FB  FR_Factor  FC  vSW  vS  vW  Re  Cf  Csf  Cy  Ay
global Var_Know  Var   Varmove  Varmove_row  Varmove_col 
global col0  Path  PathN  
% y, ydot, s(通常)是变化的全局变量
% Bs, Bv, Ba, Br, Bw, Be, Bsj, Bvj, Baj, Brj, Bwj, Bej虽然主程序中没有涉及,  
% 但它们是重要的动力响应输出变量   
global y  ydot  s   
% E不定义为全局变量,因为ode45中对E有不同的定义,容易引起误解  
% PV0(j,i,k)为长方形刚体Bk八个顶点的xyz坐标(连体坐标系中)
% PV(:,i,k,nt)为长方形刚体Bk的八个顶点在nt时刻的xyz坐标(在静止坐标系中)    
global SJK  SOK  PV0  PV   
% 运动响应的相关输出变量
global Bs  Bv  Ba  Br  Bw  Be  Angle  Bsj  Bvj  Baj  Brj  Bwj  Bej  
% 能量测试原理及能量守恒验证的相关变量
% dK_nt为nt时的系统动能对时间的导数；E_Flyl为能量测试原理计算出的系统动能； 
% E_Kin为系统的动能；E_Gra为系统的重力势能；E0_Kin为系统初动能
global dK_nt  E_Flyl  E_Kin_k  E_Gra_k  E_Kin  E_Gra  E0_Kin 
% FR[6x(Step_n+1)]为静水恢复力矩阵； 
% FDF[3x(Step_n+1)]为重力锚动摩擦力矩阵
% FSC[3x(Step_n+1)xN]为海底对海底卧链的约束力(只有z方向的)
% E_FR为静水恢复力所积蓄的弹性势能
global mu  FR  E_FR  FDF  FF  FSC
% 基于牛顿定律、Kane方法求解接点约束力的相关变量
% T_Hin为通过Kane方法计算得到的物体所受力矩：M_After为经典力学方法计算得到的物体所受力矩
% MSC是使锚链平衡所施加的力矩
global T_Hin  M_After  MSC 
% F_Con为广义约束力,F_Hin_Kane为Kane方法计算得到的接点约束力 
global F_Con  F_Hin_Kane  Mcon_3N  Fcon_3N  Mcon_k  

%　附加质量力附加阻尼力及粘性摩擦阻力相关
global  v_t  v_C  v_Cdot  F_V  T_V  F_A  T_A  F_D  T_D   F_Vx   F_Vz  F_Ax   F_Az  F_Dx   F_Dz 
% 外部约束相关变量
global  gt_xz  gtdot_xz  Error_e_xz  lamda_xz  gt_z  gtdot_z  Error_e_z  lamda_z  Bsub_xz  Csub_xz  Bsub_z  Csub_z 
 

%% ★★★定义输入参数★★★
%% 输入时间步长Delta_t、积分步数Step_n、重力加速度g、海床阻尼系数c_Seabed(暂粗略考虑)（N/m^3）-软土
% 积分时间步长会影响数值计算的收敛性                                               ============================== 弹性地基梁模型有待改进 =============================
Delta_t=1.0e-3;   Step_n=60/Delta_t;    c_Seabed=5.0e6;       
% g为重力加速度，注意其方向，因为此坐标系内z0轴竖直向上；rou_W为海水密度(kg/m^3)；
g=-9.81;                                              rou_W=1.03e3;
% 海水阻力系数C_D、海水惯性系数C_M、海水粘性系数C_F(暂粗略考虑)
C_D=2;                       C_M=2.0;        C_F=0.1; 
% ks为船舶表面粗糙度特征高度(m)(新船取150e-6)；upsilon为运动粘度；
ks=150e-6;                                   upsilon=1.0e-6; 
% 船舶停止标志位
BvSFlag=0;

%% ★多体系统刚体数N、自由度数量Free及低序体连接阵列Low(行向量)★
% N为多体系统的刚体数量； Free为系统总自由度数量；  Low为低序体连接阵列
N=72;                   Free=6*N;                Low=0:N-1;  

%% ★初始时刻多体系统相关参数的输入及悬链线形状的求解★ 
% 输入各参数的初值
% L_Cha为锚链长度(m)； h_Cha(zL)为锚链上端到海底的距离(m)；
L_Cha=15;             h_Cha=0;                          
% N_Cha为锚链划分的段数；m_CU为锚链的单位质量(kg/m)；V_CU为锚链的单位体积(m^3)；
N_Cha=5;               m_CU=143.5120;             V_CU=m_CU/7.85e3;                      
%  FB_CU为锚链单位长度的浮力；w_CUW为锚链单位长度的湿重； 
FB_CU=-rou_W*V_CU*g;       w_CUW=-m_CU*g-FB_CU;      
% L_Pre为拦阻索长度(m)； h_Pre(zL)为拦阻索上端到海底的距离(m)；x_Pre为在海底的水平投影(m)；
L_Pre=192;              h_Pre=6;                            x_Pre=187;         
% N_Pre为拦阻索划分的段数；m_PU为拦阻索单位质量(kg/m)；V_PU为拦阻索单位体积(m^3)；
N_Pre=N-3-N_Cha;         m_PU=6.9*2;                V_PU=m_PU/7.85e3;                      
% FB_PU为拦阻索单位长度的浮力；w_PUW为锚链单位长度的湿重； 
FB_PU=-rou_W*V_PU*g;       w_PUW=-m_PU*g-FB_PU;      
% 确定锚链各链环坐标及各物体初始偏角
theta(1:N_Cha+1)=0;   
X_HoC=[0 : L_Cha/N_Cha : L_Cha].'; 
Z_HoC=zeros(length(X_HoC),1);
Y_HoC=Z_HoC; 
% 确定拦阻索各点坐标及各物体初始偏角
if L_Pre^2<=x_Pre^2+h_Pre^2
    error('Wrong! The prevention cable is too short!');
elseif L_Pre>x_Pre+h_Pre
    error('Wrong! The prevention cable is too long!');
else
    % 先假定存在卧链
    CatenaryShape_HorizontalPC;
    if x_Pre-x_PreS<=1.0e-6
        CatenaryShape_NoHorizontalPC;
    end
end
% theta为多体系统各物体的初始偏角(相对于x0轴的夹角,在x0o0z0平面内)
theta(N_Cha+2:N-2)=theta_Pre;  theta(N-1)=pi/2;  theta(N)=0;

%% ★各刚体在局部连体坐标系下的几何、物理特性参数(先在连体坐标系中表示)★
% L_ChaE为各级锚链杆的长度；  R_Cha为锚链杆的等效半径
L_ChaE=L_Cha/N_Cha;      R_Cha=sqrt(m_CU/7.85e3/pi);
% L_PreU为各级缆索杆的长度；  R_Pre为缆索杆的等效半径
L_PreE=L_Pre/N_Pre;      R_Pre=sqrt(m_PU/7.85e3/pi);
% L_Aux为辅助杆件的长度; M_Aux为辅助杆件的质量（辅助杆件的半径为0.1m）
L_Aux=1;                M_Aux=10;            R_Aux=R_Pre;  
% Len1和LenN分别为重力锚和船舶的几何尺寸（长宽高）
Len1(1:3,1)=[1.4; 5.6; 2.7];  
Len_Cha=[L_ChaE; sqrt(2)*R_Cha; sqrt(2)*R_Cha];  
Len_Cha=repmat(Len_Cha,1,N_Cha);
Len_Pre=[L_PreE; sqrt(2)*R_Pre; sqrt(2)*R_Pre];  
Len_Pre=repmat(Len_Pre,1,N_Pre);
Len_Aux=[L_Aux; sqrt(2)*R_Aux; sqrt(2)*R_Aux];
LenN(1:3,1)=[40; 10; 5];
Len=[Len1  Len_Cha  Len_Pre  Len_Aux  LenN];
Rad=[0, R_Cha*ones(1,N_Cha), R_Pre*ones(1,N_Pre),R_Aux, 0]; 
% 重力锚体积
Vol_Anc=prod(Len1);
% M_Shi为船舶质量850t；        M_Anc为重力锚的质量50t；    
M_Shi=8.5e5;                  M_Anc=5.0e4;
% Mass(K)<1*N>为系统的质量矩阵(均为干重)(行向量)  
Mass(1)=M_Anc;  Mass(2:N_Cha+1)=m_CU*L_ChaE;  Mass(N_Cha+2:N-2)=m_PU*L_PreE;  Mass(N-1)=M_Aux;  Mass(N)=M_Shi;  
% Dis为船舶排水质量(kg)
Dis=-Tv_Pre/g+Mass(N-1)+Mass(N); 
% Draft为初始时刻船舶正浮状态的吃水(正数)
Draft=Dis/rou_W/prod(LenN([1 2]));  
% S_Wet为船舶湿表面
S_Wet=2*Draft*(LenN(1)+LenN(2))+prod(LenN(1:2));
% IT_Local(I,J,K)<3*3*N>为惯性并矢式阵
IT_Local(:,:,1)=diag([Len1(2)^2+Len1(3)^2,Len1(1)^2+Len1(3)^2,Len1(1)^2+Len1(2)^2])/12*Mass(1);
for k=2:N_Cha+1
    IT_Local(:,:,k)=diag([R_Cha^2/2,  L_ChaE^2/12,  L_ChaE^2/12])*Mass(k);
end
for k=N_Cha+2:N-2
    IT_Local(:,:,k)=diag([R_Pre^2/2,  L_PreE^2/12,  L_PreE^2/12])*Mass(k);
end
IT_Local(:,:,N-1)=diag([R_Aux^2/2,  L_Aux^2/12,  L_Aux^2/12])*Mass(N-1);
IT_Local(:,:,N)=diag([LenN(2)^2+LenN(3)^2,LenN(1)^2+LenN(3)^2,LenN(1)^2+LenN(2)^2])/12*Mass(N);

%% ★说明各刚体连接关系的矢量★  
% 将o_Anc坐标系变换到o0x0y0z0坐标系!☆☆☆
% Dis_A2SC为初始时刻辅助杆件末端点到惯性系原点的距离(有正负,正表示位于船舶的前部)  
Dis_A2SC=-10;  
X_HoC=X_HoC-L_Cha-x_Pre+Dis_A2SC; 
Z_HoC=Z_HoC-h_Pre-L_Aux; 
X_HoP=X_HoP-x_Pre+Dis_A2SC; 
Z_HoP=Z_HoP-h_Pre-L_Aux; 
% Q_Local(I,K)<3*N>为局部坐标系下的参考点位置矢量阵(按列向量记录各个刚体)
% Q_Local(I,K)同时也是各物体连体坐标系原点间的相对位置，
% 注意走锚时Q_Local(:,1)也不会变化☆☆☆
Q_Local(1:3,1:N)=0; 
% 船舶BN的连体坐标系与静止坐标系不重合，暂假定其质心与静止坐标系原点重合
Q_Local(1,1)=-(Len1(1)+L_Cha+x_Pre-Dis_A2SC);  Q_Local(3,1)=-h_Pre-L_Aux-R_Cha; 
Q_Local(1,2)=Len1(1);  Q_Local(3,2)=R_Cha; 
Q_Local(1,3:N_Cha+2)=L_ChaE;  
Q_Local(1,N_Cha+3:N-1)=L_PreE;   
% ☆☆☆
Q_Local(1,N)=L_Aux; 
% R_Local(I,K)<3*N>为局部坐标系下的质心位置矢量阵(按列向量记录各个刚体)  
R_Local(1:3,1:N)=0;   % 假定船舶的中心在z0=0处  
R_Local(1,1)=Len1(1)/2;  R_Local(3,1)=Len1(3)/2;
R_Local(1,2:N_Cha+1)=L_ChaE/2; 
R_Local(1,N_Cha+2:N-2)=L_PreE/2; 
R_Local(1,N-1)=L_Aux/2;
% ☆☆☆注意船舶连体坐标系的方向
R_Local(1,N)=-Dis_A2SC;  

%% ★局部坐标系及惯性坐标系下各刚体的主矢F_Local、F和主矩M_Local、M★
% 重力、浮力、水流力、海底摩擦力及海底约束力
% 多体系统中通过各刚体质心的主矢F(I,K)<3*N>和主矩M(I,K)<3*N>
% 均按列向量记录各个刚体
F(1:3,1:N)=0;        M(1:3,1:N)=0;
F_Local(1:3,1:N)=0;  M_Local(1:3,1:N)=0;
F(3,1:N)=Mass*g;			       % 所有物体的重力都必须考虑,注意此时的z0轴方向
% FB为各物体的浮力(沿z0轴向上)
FB(1:3,1:N)=0;
FB(3,1)=-rou_W*Vol_Anc*g; 
FB(3,2:N_Cha+1)=L_ChaE*FB_CU; 
FB(3,N_Cha+2:N-2)=L_PreE*FB_PU; 
FB(3,N-1)=-rou_W*pi*R_Aux^2*L_Aux*g;
FB(3,N)=-Dis*g;

% =========================================================================   ==================== 待改进 ========================== 
% 船舶的静水恢复力系数，为负值，因为吃水变化及转角变化的正负☆☆☆
zB=Draft/2;     zG=Draft;        V_Dis=Dis/rou_W; 
Ix=LenN(1)*LenN(2)^3/12;         Iy=LenN(1)^3*LenN(2)/12;
GMx=zB-zG+Ix/V_Dis;              GMy=zB-zG+Iy/V_Dis; 			
FR_Factor=zeros(6,1);
FR_Factor(3)=rou_W*g*prod(LenN([1 2]));
FR_Factor(4)=Dis*g*GMx;
FR_Factor(5)=Dis*g*GMy; 

% =========================================================================  ==================== 待改进 ========================== 
% FC为水流力,注意其正负，即方向☆☆☆
% vS为船舶的绝对速度；vW为流速(相对惯性系)；vSW为船舶相对水流的速度；
FC(1:3,1:Step_n+1,1:N)=0;  
vS(1:3,1:Step_n+1)=0;  vW(1:3,1:Step_n+1)=0;  
vS(1,1)=4.327;         vW(1,1)=1.261;   
vSW(1:3,1:Step_n+1)=0; 
vSW(1,1)=vS(1,1)-vW(1,1);    
Re=abs(vSW(1,1))*LenN(1)/upsilon;
Cf=0.075/(log(Re)-2)^2;
Csf=(105*(ks/LenN(1))^(1/3)-0.64)*1.0e-3; 
FC(1,1,N)=-sign(vSW(1,1))*(Cf+Csf)*rou_W*S_Wet*vSW(1,1).^2/2; 
% 2维时用不到                                                                 ==================== 待改进 ========================== 
Cy=0.5;  Ay=80;  
FC(2,1,N)=-sign(vSW(2,1))*rou_W*Cy*Ay*vSW(2,1).^2/2; 

% mu为海底摩擦系数
mu=0.3;  
FDF(1:3,1:Step_n+1)=0; 
FF(1:3,1:Step_n+1)=0;
% FSC为海底对各卧链的约束力(仅有z方向的)-弹性地基梁作用
FSC(1:3,1:Step_n+1,1:N)=0;  
FSC(3,1,1)=-F(3,1)-FB(3,1);
for k=2:N-2
    if theta(k)==0
        if k<=N_Cha+1
            FSC(3,1,k)=L_ChaE*w_CUW;
        else
            FSC(3,1,k)=L_PreE*w_PUW;
        end
    end
end
MSC(1:3,1:Step_n+1,1:N)=0;  



% 锚链和拦阻索的水流力
F_V=zeros(Step_n+1,N);
T_V=zeros(Step_n+1,N);
F_A=zeros(Step_n+1,N);
T_A=zeros(Step_n+1,N);
F_D=zeros(Step_n+1,N);
T_D=zeros(Step_n+1,N);
F_Vx=zeros(Step_n+1,N);
F_Vz=zeros(Step_n+1,N);
F_Ax=zeros(Step_n+1,N);
F_Az=zeros(Step_n+1,N);
F_Dx=zeros(Step_n+1,N);
F_Dz=zeros(Step_n+1,N);


% =========================================================================  ====================== 待改进 ========================== 
T_Hin(1:3,1:Step_n+1,1:N)=0; 
% F_Hin_Kane为接点约束力
F_Hin_Kane(1:3,1:Step_n+1,N)=0; 

%% ★记入联接点的物体特性★
% 关于接点约束及广义坐标、广义速度等的定义
% Var_Know表示已被约束固定的自由度编号，无需求解(是一行向量)
% Var为待求解自由度的编号(是一行向量)
% 缩减自由度变量得到须积分求解的自由度
Var_Know=1:6*N;  
% 不考虑重力锚的转动，只考虑其x方向的移动
Var=[4:3*N, 3*N+1:3*N+3];  
% Var=[5:3:3*N, 3*N+1]; 
Var_Know(Var)=[];
% ☆☆☆注意
% Varmove为一向量，其元素表示待求解的移动自由度编号，是移动位移而非移动速率！☆☆☆
% 单独列出相对位移矢量 s(1:3,1)对应的分量号及刚体号，待求解的自由度
Varmove=3*N+1:3*N+3;  
Varmove_row=[1 2 3];  Varmove_col=[1 1 1];   

%% ★记入各待定变量及其导数的初始值--定义初始状态！★
% 已知广义变量（速度和加速度）的初值,y表示广义速率(为一列向量)
% 其他广义坐标及其导数初始值
y(1:Free,1)=0;  
% omiga0_A2P为辅助杆件相对其内接刚体的角速度分量（绕y(N-1)轴）
vs0x=vS(1,1);  omiga0_A2P=vs0x/L_Aux;  
% 初动能 =（辅助杆件质心的平动动能+其绕质心的转动动能）+船舶的平动动能
E0_Kin=Mass(N-1)*(vs0x/2)^2/2+IT_Local(2,2,N-1)*omiga0_A2P^2/2+Mass(N)*vs0x^2/2;
% omiga0_S2A为船舶相对于其内接刚体(辅助竖直杆件)的角速度
omiga0_S2A=-omiga0_A2P;
y(3*(N-1)-1)=omiga0_A2P;  
y(3*N-1)=omiga0_S2A;  
% 广义速率的导数就是加速度(为一列向量)
ydot(1:Free,1)=0;	
% 广义速率的相对速度部分(按列向量记录各个刚体)
s(1:3,1:N)=0;
y0=y;  ydot0=ydot;  s0=s;

%% 由初始基矢量位置求初始变换矩阵
% n0(:,:,k)=[e1k  e2k  e3k].'=[e1k.';e2k.';e3k.']为刚体Bk连体基矢量组成的矩阵
% n0(:,:,k).'是刚体Bk的连体坐标系向静止惯性系转换的变换矩阵
for k=1:N
    n0(:,:,k)=[cos(theta(k)),  0,  sin(theta(k));
                           0,  1,              0;
               -sin(theta(k)), 0,  cos(theta(k))];
end
% SJK0表示初始时刻K系向J系(相邻刚体间)的变换矩阵 
% SJK0(:,:,k)表示初始时刻刚体Bk向其内接刚体Bj转换的变换矩阵
E=eye(3,3);
SJK0(:,:,1)=E*n0(:,:,1)';
for k=2:N
    SJK0(:,:,k)=n0(:,:,Low(k))*n0(:,:,k)';
end

%% 欧拉参数及其导数的初始值(按列向量记录各个刚体间的变换)
% epsi0为初始时刻的欧拉参数
for k=1:N
   epsi0(4,k)=sqrt(1+trace(SJK0(:,:,k)))/2;
   epsi0(1:3,k)=real((sqrt((1+diag(SJK0(:,:,k)))/2-epsi0(4,k)^2)));
end
% 怎么决定正负？回带！  尤其要小心正负号，出错不止一次
for k=1:N-1
    epsi0(2,k)=-epsi0(2,k);
end
for k=1:N
      for i=1:3
         for j=1:3
             SJK0_epsi0temp(i,j)=epsi0(i,k)*epsi0(j,k);
         end
      end
      SJK0_epsi0(:,:,k)=2*(SJK0_epsi0temp+epsi0(4,k)*(epsi0(4,k)*E+Dual_Mat(epsi0(1:3,k))))-E;
end 
% 回带验证欧拉参数的正负号
for k=1:N
      for i=1:3
         for j=1:3
             if abs(SJK0_epsi0(i,j,k)-SJK0(i,j,k))>1.0e-6
                 error([ 'There are something wrong on the symbols of Euler Parameters on the', num2str(k),'th body !'])
             end
         end
      end
end 
% 积分变量初值的列向量
temp=[];
for i=1:length(Varmove)   
    temp=[temp;s0(Varmove_row(i),Varmove_col(i))];	% 由此可看出两个变量的意义
end
col0=[y0(Var);reshape(epsi0,4*N,1);temp];
%% ★★★输入参数读取结束，下面开始相关计算★★★

%% 通路体阵列
% PathN行向量表示了各刚体进行Low算子计算通向静止坐标系B0的次数
Path(1:N,1:N)=0;
PathN(1:N)=N;   %初始化（不可置零）
Path(1,:)=1:N;
for k=1:N
    for i=2:N
        if Low(Path(i-1,k))==0
            PathN(k)=i-1;
            break
        end
        Path(i,k)=Low(Path(i-1,k));
    end
end

%% 初始顶点位置（仅对长方体）
for k=1:N
    for i=1:8
        for j=1:3
            PV0(j,i,k)=Len(j,k)/2*(-1)^(fix((i-1)/(2^(j-1))));
        end
    end
end
num=[1 2 6 8 7 5 6 2 3 7 6 2 4 8 5 1 3 4];	%用于画长方体的顶点序列，水平面有交叉斜线

%% ★求解一阶微分方程★
tic;  disp(' ');  
disp('    Kane equations are calculating ! ');  
% 下面语句等效于[T,Y]=solver('odefun',tspan,y0,options); 
[T,Y]=Initial_Fixed_ode45('Dyn_Res_of_Anc_Cha_Anc1_Fun');  
disp(['    Time for calculating is ',num2str(toc),' s ']);  
time_c=toc;

t=Delta_t*(0:size(Bs,2)-1);
% 动摩擦力作负功(FDF始终为负值)；水流力做功可正可负(FC时正时负) 
for nt=1:length(t) 
    dW_FF(nt)=FDF(1,nt)*Bv(1,nt,1)*Delta_t;
    dW_FC(nt)=FC(1,nt,N)*Bv(1,nt,N)*Delta_t;
    W_FF(nt)=sum(dW_FF(1:nt));
    W_FC(nt)=sum(dW_FC(1:nt));

    for k=1:N
        if k==1
            W_FSCk(nt,k)=FSC(3,nt,k)*(Bs(3,nt,k)+h_Pre+L_Aux+R_Cha-Len1(3)/2)/2;
        else
            W_FSCk(nt,k)=FSC(3,nt,k)*(Bs(3,nt,k)+h_Pre+L_Aux)/2;
        end
        
        
        dW_FVk(nt,k)=F_Vx(nt,k)*Bv(1,nt,k)*Delta_t+F_Vz(nt,k)*Bv(3,nt,k)*Delta_t;
        dW_FAk(nt,k)=F_Ax(nt,k)*Bv(1,nt,k)*Delta_t+F_Az(nt,k)*Bv(3,nt,k)*Delta_t;
        dW_FDk(nt,k)=F_Dx(nt,k)*Bv(1,nt,k)*Delta_t+F_Dz(nt,k)*Bv(3,nt,k)*Delta_t;
        dW_TAk(nt,k)=T_A(nt,k)*Bw(2,nt,k)*Delta_t;
        dW_TDk(nt,k)=T_D(nt,k)*Bw(2,nt,k)*Delta_t;
        
        
    end
    W_FSC(nt)=sum(W_FSCk(nt,1:N));  
    
    % 注意：与FSC不同。dW_FVk(nt,k)表示Delta_t时间内第k个物体的切向水流力所做的功；
    % dW_FV(nt)表示Delta_t时间内切向力水流力对系统所做的功；
    % W_FV(nt)表示nt时刻切向水流力对系统所做的功。
    dW_FV(nt)=sum(dW_FVk(nt,1:N));   W_FV(nt)=sum(dW_FV(1:nt)); 
    dW_FA(nt)=sum(dW_FAk(nt,1:N));   W_FA(nt)=sum(dW_FA(1:nt));  
    dW_FD(nt)=sum(dW_FDk(nt,1:N));   W_FD(nt)=sum(dW_FD(1:nt));
    dW_TA(nt)=sum(dW_TAk(nt,1:N));   W_TA(nt)=sum(dW_TA(1:nt));
    dW_TD(nt)=sum(dW_TDk(nt,1:N));   W_TD(nt)=sum(dW_TD(1:nt));  
end  

E_all=E_Kin+E_Gra+E_FR-W_FF-W_FC-W_FSC-W_FV-W_FA-W_FD-W_TA-W_TD; 

%% 存储数据
save 1.mat

%% 动画演示
set(0,'DefaultFigureColor','w');
set(0,'DefaultAxesFontname','Times New Roman');
color=['k','g','b','r','m'];  color=repmat(color,1,N);
fig1=figure(1);  view([30 30]);  hold on;  zoom;
axis('equal'); 
axis([-Len1(1)-L_Cha-x_Pre-LenN(1) 100, -100 100, -h_Pre-R_Cha  2*h_Pre]);
set(fig1,'Position',[100 65 800 600]);  
for i=1:N
    hand(i)=line('color',color(i),'linewidth',2,'linestyle','-','erasemode','xor');
end
% 每1s画一次
% for nt=1:max(1,0.1/Delta_t):Step_n
for nt=1: Step_n/1.0e2:length(t)
    for k=1:N
        set(hand(k),'xdata',PV(1,num,k,nt),'ydata',PV(2,num,k,nt),'zdata',PV(3,num,k,nt));
    end
    drawnow;
    pause(0.1);
end
xlabel('x');  ylabel('y');  zlabel('z');   

%% 画图-能量测试原理验证和能量守恒验证
figure;  plot(t,E_Kin,'or',t,E_Flyl,'.b','Markersize',2);
title('Testing the accuracy of numerical simulations');  box off;
xlabel('Time (s)');  ylabel('Kinetic Energy (J)');  

figure;  plot(t,E_all,'-k'); 
title('Testing based on conservation principle of energy');  box off;
xlabel('Time (s)');  ylabel('Total Energy (J)');

figure; plot(t,Bv(1,1:length(t),1),'k',t,Bv(1,1:length(t),N),'b','Linewidth',2);

figure; plot(t,sqrt(F_Hin_Kane(1,1:length(t),2).^2+F_Hin_Kane(3,1:length(t),2).^2),'r')



 