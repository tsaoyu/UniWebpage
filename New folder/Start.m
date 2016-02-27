
% ȫ�׳����ǻ���Kane�������ê���Ķ�����Ӧ(2ά���տ�ʼ����ê������Ȼ����ê)
% �����������Ҫ�����Ǹ�������������,������ֺ����ĳ�ֵ
% Ȼ����������������һ�׳�΢�ַ����飬�����ж�����ʾ

% ע�⣺�ڴ�����ê������֮�����һ����ֱ�˼��Ա㶨�崬���ĳ��ٶȣ�����
% ����êê�鱻��Ϊ�����壬����߷ֱ�Ϊ1.4m*5.6m*2.7m������12.5tΪ1.4m*1.4m*2.7m��
% ê��Ϊ�޵�ê������������Ϊ60mm,����ê���ĵ�λ����m_CUΪ2*71.756=143.5120kg/m, ʪ��w_CUW=-m_CU*g-FB_CU  
% ê��������������ΪԲ����˼�����Ч�뾶R_CE=sqrt(m_CU/7.85e3/pi)=0.0763m  
% ��������ЧΪ��˿������˿ֱ��7mm��19������ģ���еĵ�Ч�뾶R_PE=sqrt(m_PU/7.85e3/pi)=0.0167�� (2) 0.0035m*sqrt(38)=0.0216m
% ��������Ϊ�����壬����߷ֱ�Ϊ40m*10m*5m��500tʱ��ˮ1.25m,850tʱ��ˮ2.125m��

% ������ֹ����ϵ���̶���o0Ϊԭ�㣬x0��ˮƽ���ң�z0����ֱ���ϣ�y0����ֱ����
% ��ֹ����ϵԭ��λ�ڴ�������״̬��ˮ�������Ĵ����ɲ��������غϣ�
% �Ժ�����Ϊ�������ܵ���������

clc;  clear ;  close all;

%% ����ȫ�ֱ���
global Delta_t  Step_n  c_Seabed  g  rou_W  C_D  C_M  C_F  ks  upsilon  BvSFlag 
global N  Free  Low
global N_Cha  L_Pre  h_Pre  x_Pre  w_PUW  N_Pre  X_HoP  Z_HoP  Tv_Pre  theta  
global R_Cha  Len1  Len_Cha  Len_Pre  L_Aux  LenN  Len  Rad  Vol_Anc  M_Shi  M_Anc  Mass 
global Dis  S_Wet  IT_Local Q_Local  R_Local      
global F  M  F_Local  M_Local  FB  FR_Factor  FC  vSW  vS  vW  Re  Cf  Csf  Cy  Ay
global Var_Know  Var   Varmove  Varmove_row  Varmove_col 
global col0  Path  PathN  
% y, ydot, s(ͨ��)�Ǳ仯��ȫ�ֱ���
% Bs, Bv, Ba, Br, Bw, Be, Bsj, Bvj, Baj, Brj, Bwj, Bej��Ȼ��������û���漰,  
% ����������Ҫ�Ķ�����Ӧ�������   
global y  ydot  s   
% E������Ϊȫ�ֱ���,��Ϊode45�ж�E�в�ͬ�Ķ���,�����������  
% PV0(j,i,k)Ϊ�����θ���Bk�˸������xyz����(��������ϵ��)
% PV(:,i,k,nt)Ϊ�����θ���Bk�İ˸�������ntʱ�̵�xyz����(�ھ�ֹ����ϵ��)    
global SJK  SOK  PV0  PV   
% �˶���Ӧ������������
global Bs  Bv  Ba  Br  Bw  Be  Angle  Bsj  Bvj  Baj  Brj  Bwj  Bej  
% ��������ԭ�������غ���֤����ر���
% dK_ntΪntʱ��ϵͳ���ܶ�ʱ��ĵ�����E_FlylΪ��������ԭ��������ϵͳ���ܣ� 
% E_KinΪϵͳ�Ķ��ܣ�E_GraΪϵͳ���������ܣ�E0_KinΪϵͳ������
global dK_nt  E_Flyl  E_Kin_k  E_Gra_k  E_Kin  E_Gra  E0_Kin 
% FR[6x(Step_n+1)]Ϊ��ˮ�ָ������� 
% FDF[3x(Step_n+1)]Ϊ����ê��Ħ��������
% FSC[3x(Step_n+1)xN]Ϊ���׶Ժ���������Լ����(ֻ��z�����)
% E_FRΪ��ˮ�ָ���������ĵ�������
global mu  FR  E_FR  FDF  FF  FSC
% ����ţ�ٶ��ɡ�Kane�������ӵ�Լ��������ر���
% T_HinΪͨ��Kane��������õ��������������أ�M_AfterΪ������ѧ��������õ���������������
% MSC��ʹê��ƽ����ʩ�ӵ�����
global T_Hin  M_After  MSC 
% F_ConΪ����Լ����,F_Hin_KaneΪKane��������õ��Ľӵ�Լ���� 
global F_Con  F_Hin_Kane  Mcon_3N  Fcon_3N  Mcon_k  

%������������������������ճ��Ħ���������
global  v_t  v_C  v_Cdot  F_V  T_V  F_A  T_A  F_D  T_D   F_Vx   F_Vz  F_Ax   F_Az  F_Dx   F_Dz 
% �ⲿԼ����ر���
global  gt_xz  gtdot_xz  Error_e_xz  lamda_xz  gt_z  gtdot_z  Error_e_z  lamda_z  Bsub_xz  Csub_xz  Bsub_z  Csub_z 
 

%% ���ﶨ�������������
%% ����ʱ�䲽��Delta_t�����ֲ���Step_n���������ٶ�g����������ϵ��c_Seabed(�ݴ��Կ���)��N/m^3��-����
% ����ʱ�䲽����Ӱ����ֵ�����������                                               ============================== ���Եػ���ģ���д��Ľ� =============================
Delta_t=1.0e-3;   Step_n=60/Delta_t;    c_Seabed=5.0e6;       
% gΪ�������ٶȣ�ע���䷽����Ϊ������ϵ��z0����ֱ���ϣ�rou_WΪ��ˮ�ܶ�(kg/m^3)��
g=-9.81;                                              rou_W=1.03e3;
% ��ˮ����ϵ��C_D����ˮ����ϵ��C_M����ˮճ��ϵ��C_F(�ݴ��Կ���)
C_D=2;                       C_M=2.0;        C_F=0.1; 
% ksΪ��������ֲڶ������߶�(m)(�´�ȡ150e-6)��upsilonΪ�˶�ճ�ȣ�
ks=150e-6;                                   upsilon=1.0e-6; 
% ����ֹͣ��־λ
BvSFlag=0;

%% �����ϵͳ������N�����ɶ�����Free����������������Low(������)��
% NΪ����ϵͳ�ĸ��������� FreeΪϵͳ�����ɶ�������  LowΪ��������������
N=72;                   Free=6*N;                Low=0:N-1;  

%% ���ʼʱ�̶���ϵͳ��ز��������뼰��������״������ 
% ����������ĳ�ֵ
% L_ChaΪê������(m)�� h_Cha(zL)Ϊê���϶˵����׵ľ���(m)��
L_Cha=15;             h_Cha=0;                          
% N_ChaΪê�����ֵĶ�����m_CUΪê���ĵ�λ����(kg/m)��V_CUΪê���ĵ�λ���(m^3)��
N_Cha=5;               m_CU=143.5120;             V_CU=m_CU/7.85e3;                      
%  FB_CUΪê����λ���ȵĸ�����w_CUWΪê����λ���ȵ�ʪ�أ� 
FB_CU=-rou_W*V_CU*g;       w_CUW=-m_CU*g-FB_CU;      
% L_PreΪ����������(m)�� h_Pre(zL)Ϊ�������϶˵����׵ľ���(m)��x_PreΪ�ں��׵�ˮƽͶӰ(m)��
L_Pre=192;              h_Pre=6;                            x_Pre=187;         
% N_PreΪ���������ֵĶ�����m_PUΪ��������λ����(kg/m)��V_PUΪ��������λ���(m^3)��
N_Pre=N-3-N_Cha;         m_PU=6.9*2;                V_PU=m_PU/7.85e3;                      
% FB_PUΪ��������λ���ȵĸ�����w_PUWΪê����λ���ȵ�ʪ�أ� 
FB_PU=-rou_W*V_PU*g;       w_PUW=-m_PU*g-FB_PU;      
% ȷ��ê�����������꼰�������ʼƫ��
theta(1:N_Cha+1)=0;   
X_HoC=[0 : L_Cha/N_Cha : L_Cha].'; 
Z_HoC=zeros(length(X_HoC),1);
Y_HoC=Z_HoC; 
% ȷ���������������꼰�������ʼƫ��
if L_Pre^2<=x_Pre^2+h_Pre^2
    error('Wrong! The prevention cable is too short!');
elseif L_Pre>x_Pre+h_Pre
    error('Wrong! The prevention cable is too long!');
else
    % �ȼٶ���������
    CatenaryShape_HorizontalPC;
    if x_Pre-x_PreS<=1.0e-6
        CatenaryShape_NoHorizontalPC;
    end
end
% thetaΪ����ϵͳ������ĳ�ʼƫ��(�����x0��ļн�,��x0o0z0ƽ����)
theta(N_Cha+2:N-2)=theta_Pre;  theta(N-1)=pi/2;  theta(N)=0;

%% ��������ھֲ���������ϵ�µļ��Ρ��������Բ���(������������ϵ�б�ʾ)��
% L_ChaEΪ����ê���˵ĳ��ȣ�  R_ChaΪê���˵ĵ�Ч�뾶
L_ChaE=L_Cha/N_Cha;      R_Cha=sqrt(m_CU/7.85e3/pi);
% L_PreUΪ���������˵ĳ��ȣ�  R_PreΪ�����˵ĵ�Ч�뾶
L_PreE=L_Pre/N_Pre;      R_Pre=sqrt(m_PU/7.85e3/pi);
% L_AuxΪ�����˼��ĳ���; M_AuxΪ�����˼��������������˼��İ뾶Ϊ0.1m��
L_Aux=1;                M_Aux=10;            R_Aux=R_Pre;  
% Len1��LenN�ֱ�Ϊ����ê�ʹ����ļ��γߴ磨����ߣ�
Len1(1:3,1)=[1.4; 5.6; 2.7];  
Len_Cha=[L_ChaE; sqrt(2)*R_Cha; sqrt(2)*R_Cha];  
Len_Cha=repmat(Len_Cha,1,N_Cha);
Len_Pre=[L_PreE; sqrt(2)*R_Pre; sqrt(2)*R_Pre];  
Len_Pre=repmat(Len_Pre,1,N_Pre);
Len_Aux=[L_Aux; sqrt(2)*R_Aux; sqrt(2)*R_Aux];
LenN(1:3,1)=[40; 10; 5];
Len=[Len1  Len_Cha  Len_Pre  Len_Aux  LenN];
Rad=[0, R_Cha*ones(1,N_Cha), R_Pre*ones(1,N_Pre),R_Aux, 0]; 
% ����ê���
Vol_Anc=prod(Len1);
% M_ShiΪ��������850t��        M_AncΪ����ê������50t��    
M_Shi=8.5e5;                  M_Anc=5.0e4;
% Mass(K)<1*N>Ϊϵͳ����������(��Ϊ����)(������)  
Mass(1)=M_Anc;  Mass(2:N_Cha+1)=m_CU*L_ChaE;  Mass(N_Cha+2:N-2)=m_PU*L_PreE;  Mass(N-1)=M_Aux;  Mass(N)=M_Shi;  
% DisΪ������ˮ����(kg)
Dis=-Tv_Pre/g+Mass(N-1)+Mass(N); 
% DraftΪ��ʼʱ�̴�������״̬�ĳ�ˮ(����)
Draft=Dis/rou_W/prod(LenN([1 2]));  
% S_WetΪ����ʪ����
S_Wet=2*Draft*(LenN(1)+LenN(2))+prod(LenN(1:2));
% IT_Local(I,J,K)<3*3*N>Ϊ���Բ�ʸʽ��
IT_Local(:,:,1)=diag([Len1(2)^2+Len1(3)^2,Len1(1)^2+Len1(3)^2,Len1(1)^2+Len1(2)^2])/12*Mass(1);
for k=2:N_Cha+1
    IT_Local(:,:,k)=diag([R_Cha^2/2,  L_ChaE^2/12,  L_ChaE^2/12])*Mass(k);
end
for k=N_Cha+2:N-2
    IT_Local(:,:,k)=diag([R_Pre^2/2,  L_PreE^2/12,  L_PreE^2/12])*Mass(k);
end
IT_Local(:,:,N-1)=diag([R_Aux^2/2,  L_Aux^2/12,  L_Aux^2/12])*Mass(N-1);
IT_Local(:,:,N)=diag([LenN(2)^2+LenN(3)^2,LenN(1)^2+LenN(3)^2,LenN(1)^2+LenN(2)^2])/12*Mass(N);

%% ��˵�����������ӹ�ϵ��ʸ����  
% ��o_Anc����ϵ�任��o0x0y0z0����ϵ!����
% Dis_A2SCΪ��ʼʱ�̸����˼�ĩ�˵㵽����ϵԭ��ľ���(������,����ʾλ�ڴ�����ǰ��)  
Dis_A2SC=-10;  
X_HoC=X_HoC-L_Cha-x_Pre+Dis_A2SC; 
Z_HoC=Z_HoC-h_Pre-L_Aux; 
X_HoP=X_HoP-x_Pre+Dis_A2SC; 
Z_HoP=Z_HoP-h_Pre-L_Aux; 
% Q_Local(I,K)<3*N>Ϊ�ֲ�����ϵ�µĲο���λ��ʸ����(����������¼��������)
% Q_Local(I,K)ͬʱҲ�Ǹ�������������ϵԭ�������λ�ã�
% ע����êʱQ_Local(:,1)Ҳ����仯����
Q_Local(1:3,1:N)=0; 
% ����BN����������ϵ�뾲ֹ����ϵ���غϣ��ݼٶ��������뾲ֹ����ϵԭ���غ�
Q_Local(1,1)=-(Len1(1)+L_Cha+x_Pre-Dis_A2SC);  Q_Local(3,1)=-h_Pre-L_Aux-R_Cha; 
Q_Local(1,2)=Len1(1);  Q_Local(3,2)=R_Cha; 
Q_Local(1,3:N_Cha+2)=L_ChaE;  
Q_Local(1,N_Cha+3:N-1)=L_PreE;   
% ����
Q_Local(1,N)=L_Aux; 
% R_Local(I,K)<3*N>Ϊ�ֲ�����ϵ�µ�����λ��ʸ����(����������¼��������)  
R_Local(1:3,1:N)=0;   % �ٶ�������������z0=0��  
R_Local(1,1)=Len1(1)/2;  R_Local(3,1)=Len1(3)/2;
R_Local(1,2:N_Cha+1)=L_ChaE/2; 
R_Local(1,N_Cha+2:N-2)=L_PreE/2; 
R_Local(1,N-1)=L_Aux/2;
% ����ע�⴬����������ϵ�ķ���
R_Local(1,N)=-Dis_A2SC;  

%% ��ֲ�����ϵ����������ϵ�¸��������ʸF_Local��F������M_Local��M��
% ������������ˮ����������Ħ����������Լ����
% ����ϵͳ��ͨ�����������ĵ���ʸF(I,K)<3*N>������M(I,K)<3*N>
% ������������¼��������
F(1:3,1:N)=0;        M(1:3,1:N)=0;
F_Local(1:3,1:N)=0;  M_Local(1:3,1:N)=0;
F(3,1:N)=Mass*g;			       % ������������������뿼��,ע���ʱ��z0�᷽��
% FBΪ������ĸ���(��z0������)
FB(1:3,1:N)=0;
FB(3,1)=-rou_W*Vol_Anc*g; 
FB(3,2:N_Cha+1)=L_ChaE*FB_CU; 
FB(3,N_Cha+2:N-2)=L_PreE*FB_PU; 
FB(3,N-1)=-rou_W*pi*R_Aux^2*L_Aux*g;
FB(3,N)=-Dis*g;

% =========================================================================   ==================== ���Ľ� ========================== 
% �����ľ�ˮ�ָ���ϵ����Ϊ��ֵ����Ϊ��ˮ�仯��ת�Ǳ仯����������
zB=Draft/2;     zG=Draft;        V_Dis=Dis/rou_W; 
Ix=LenN(1)*LenN(2)^3/12;         Iy=LenN(1)^3*LenN(2)/12;
GMx=zB-zG+Ix/V_Dis;              GMy=zB-zG+Iy/V_Dis; 			
FR_Factor=zeros(6,1);
FR_Factor(3)=rou_W*g*prod(LenN([1 2]));
FR_Factor(4)=Dis*g*GMx;
FR_Factor(5)=Dis*g*GMy; 

% =========================================================================  ==================== ���Ľ� ========================== 
% FCΪˮ����,ע�������������������
% vSΪ�����ľ����ٶȣ�vWΪ����(��Թ���ϵ)��vSWΪ�������ˮ�����ٶȣ�
FC(1:3,1:Step_n+1,1:N)=0;  
vS(1:3,1:Step_n+1)=0;  vW(1:3,1:Step_n+1)=0;  
vS(1,1)=4.327;         vW(1,1)=1.261;   
vSW(1:3,1:Step_n+1)=0; 
vSW(1,1)=vS(1,1)-vW(1,1);    
Re=abs(vSW(1,1))*LenN(1)/upsilon;
Cf=0.075/(log(Re)-2)^2;
Csf=(105*(ks/LenN(1))^(1/3)-0.64)*1.0e-3; 
FC(1,1,N)=-sign(vSW(1,1))*(Cf+Csf)*rou_W*S_Wet*vSW(1,1).^2/2; 
% 2άʱ�ò���                                                                 ==================== ���Ľ� ========================== 
Cy=0.5;  Ay=80;  
FC(2,1,N)=-sign(vSW(2,1))*rou_W*Cy*Ay*vSW(2,1).^2/2; 

% muΪ����Ħ��ϵ��
mu=0.3;  
FDF(1:3,1:Step_n+1)=0; 
FF(1:3,1:Step_n+1)=0;
% FSCΪ���׶Ը�������Լ����(����z�����)-���Եػ�������
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



% ê������������ˮ����
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


% =========================================================================  ====================== ���Ľ� ========================== 
T_Hin(1:3,1:Step_n+1,1:N)=0; 
% F_Hin_KaneΪ�ӵ�Լ����
F_Hin_Kane(1:3,1:Step_n+1,N)=0; 

%% ��������ӵ���������ԡ�
% ���ڽӵ�Լ�����������ꡢ�����ٶȵȵĶ���
% Var_Know��ʾ�ѱ�Լ���̶������ɶȱ�ţ��������(��һ������)
% VarΪ��������ɶȵı��(��һ������)
% �������ɶȱ����õ�������������ɶ�
Var_Know=1:6*N;  
% ����������ê��ת����ֻ������x������ƶ�
Var=[4:3*N, 3*N+1:3*N+3];  
% Var=[5:3:3*N, 3*N+1]; 
Var_Know(Var)=[];
% ����ע��
% VarmoveΪһ��������Ԫ�ر�ʾ�������ƶ����ɶȱ�ţ����ƶ�λ�ƶ����ƶ����ʣ�����
% �����г����λ��ʸ�� s(1:3,1)��Ӧ�ķ����ż�����ţ����������ɶ�
Varmove=3*N+1:3*N+3;  
Varmove_row=[1 2 3];  Varmove_col=[1 1 1];   

%% �����������������䵼���ĳ�ʼֵ--�����ʼ״̬����
% ��֪����������ٶȺͼ��ٶȣ��ĳ�ֵ,y��ʾ��������(Ϊһ������)
% �����������꼰�䵼����ʼֵ
y(1:Free,1)=0;  
% omiga0_A2PΪ�����˼�������ڽӸ���Ľ��ٶȷ�������y(N-1)�ᣩ
vs0x=vS(1,1);  omiga0_A2P=vs0x/L_Aux;  
% ������ =�������˼����ĵ�ƽ������+�������ĵ�ת�����ܣ�+������ƽ������
E0_Kin=Mass(N-1)*(vs0x/2)^2/2+IT_Local(2,2,N-1)*omiga0_A2P^2/2+Mass(N)*vs0x^2/2;
% omiga0_S2AΪ������������ڽӸ���(������ֱ�˼�)�Ľ��ٶ�
omiga0_S2A=-omiga0_A2P;
y(3*(N-1)-1)=omiga0_A2P;  
y(3*N-1)=omiga0_S2A;  
% �������ʵĵ������Ǽ��ٶ�(Ϊһ������)
ydot(1:Free,1)=0;	
% �������ʵ�����ٶȲ���(����������¼��������)
s(1:3,1:N)=0;
y0=y;  ydot0=ydot;  s0=s;

%% �ɳ�ʼ��ʸ��λ�����ʼ�任����
% n0(:,:,k)=[e1k  e2k  e3k].'=[e1k.';e2k.';e3k.']Ϊ����Bk�����ʸ����ɵľ���
% n0(:,:,k).'�Ǹ���Bk����������ϵ��ֹ����ϵת���ı任����
for k=1:N
    n0(:,:,k)=[cos(theta(k)),  0,  sin(theta(k));
                           0,  1,              0;
               -sin(theta(k)), 0,  cos(theta(k))];
end
% SJK0��ʾ��ʼʱ��Kϵ��Jϵ(���ڸ����)�ı任���� 
% SJK0(:,:,k)��ʾ��ʼʱ�̸���Bk�����ڽӸ���Bjת���ı任����
E=eye(3,3);
SJK0(:,:,1)=E*n0(:,:,1)';
for k=2:N
    SJK0(:,:,k)=n0(:,:,Low(k))*n0(:,:,k)';
end

%% ŷ���������䵼���ĳ�ʼֵ(����������¼���������ı任)
% epsi0Ϊ��ʼʱ�̵�ŷ������
for k=1:N
   epsi0(4,k)=sqrt(1+trace(SJK0(:,:,k)))/2;
   epsi0(1:3,k)=real((sqrt((1+diag(SJK0(:,:,k)))/2-epsi0(4,k)^2)));
end
% ��ô�����������ش���  ����ҪС�������ţ�����ֹһ��
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
% �ش���֤ŷ��������������
for k=1:N
      for i=1:3
         for j=1:3
             if abs(SJK0_epsi0(i,j,k)-SJK0(i,j,k))>1.0e-6
                 error([ 'There are something wrong on the symbols of Euler Parameters on the', num2str(k),'th body !'])
             end
         end
      end
end 
% ���ֱ�����ֵ��������
temp=[];
for i=1:length(Varmove)   
    temp=[temp;s0(Varmove_row(i),Varmove_col(i))];	% �ɴ˿ɿ�����������������
end
col0=[y0(Var);reshape(epsi0,4*N,1);temp];
%% �������������ȡ���������濪ʼ��ؼ������

%% ͨ·������
% PathN��������ʾ�˸��������Low���Ӽ���ͨ��ֹ����ϵB0�Ĵ���
Path(1:N,1:N)=0;
PathN(1:N)=N;   %��ʼ�����������㣩
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

%% ��ʼ����λ�ã����Գ����壩
for k=1:N
    for i=1:8
        for j=1:3
            PV0(j,i,k)=Len(j,k)/2*(-1)^(fix((i-1)/(2^(j-1))));
        end
    end
end
num=[1 2 6 8 7 5 6 2 3 7 6 2 4 8 5 1 3 4];	%���ڻ�������Ķ������У�ˮƽ���н���б��

%% �����һ��΢�ַ��̡�
tic;  disp(' ');  
disp('    Kane equations are calculating ! ');  
% ��������Ч��[T,Y]=solver('odefun',tspan,y0,options); 
[T,Y]=Initial_Fixed_ode45('Dyn_Res_of_Anc_Cha_Anc1_Fun');  
disp(['    Time for calculating is ',num2str(toc),' s ']);  
time_c=toc;

t=Delta_t*(0:size(Bs,2)-1);
% ��Ħ����������(FDFʼ��Ϊ��ֵ)��ˮ�������������ɸ�(FCʱ��ʱ��) 
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
    
    % ע�⣺��FSC��ͬ��dW_FVk(nt,k)��ʾDelta_tʱ���ڵ�k�����������ˮ���������Ĺ���
    % dW_FV(nt)��ʾDelta_tʱ����������ˮ������ϵͳ�����Ĺ���
    % W_FV(nt)��ʾntʱ������ˮ������ϵͳ�����Ĺ���
    dW_FV(nt)=sum(dW_FVk(nt,1:N));   W_FV(nt)=sum(dW_FV(1:nt)); 
    dW_FA(nt)=sum(dW_FAk(nt,1:N));   W_FA(nt)=sum(dW_FA(1:nt));  
    dW_FD(nt)=sum(dW_FDk(nt,1:N));   W_FD(nt)=sum(dW_FD(1:nt));
    dW_TA(nt)=sum(dW_TAk(nt,1:N));   W_TA(nt)=sum(dW_TA(1:nt));
    dW_TD(nt)=sum(dW_TDk(nt,1:N));   W_TD(nt)=sum(dW_TD(1:nt));  
end  

E_all=E_Kin+E_Gra+E_FR-W_FF-W_FC-W_FSC-W_FV-W_FA-W_FD-W_TA-W_TD; 

%% �洢����
save 1.mat

%% ������ʾ
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
% ÿ1s��һ��
% for nt=1:max(1,0.1/Delta_t):Step_n
for nt=1: Step_n/1.0e2:length(t)
    for k=1:N
        set(hand(k),'xdata',PV(1,num,k,nt),'ydata',PV(2,num,k,nt),'zdata',PV(3,num,k,nt));
    end
    drawnow;
    pause(0.1);
end
xlabel('x');  ylabel('y');  zlabel('z');   

%% ��ͼ-��������ԭ����֤�������غ���֤
figure;  plot(t,E_Kin,'or',t,E_Flyl,'.b','Markersize',2);
title('Testing the accuracy of numerical simulations');  box off;
xlabel('Time (s)');  ylabel('Kinetic Energy (J)');  

figure;  plot(t,E_all,'-k'); 
title('Testing based on conservation principle of energy');  box off;
xlabel('Time (s)');  ylabel('Total Energy (J)');

figure; plot(t,Bv(1,1:length(t),1),'k',t,Bv(1,1:length(t),N),'b','Linewidth',2);

figure; plot(t,sqrt(F_Hin_Kane(1,1:length(t),2).^2+F_Hin_Kane(3,1:length(t),2).^2),'r')



 