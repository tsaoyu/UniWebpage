function [out1,out2,out3] =Dyn_Res_of_Anc_Cha_Anc1_Fun(t,col,flag)
% 本函数为锚链多体系统的龙格—库塔积分函数(缩减的基本运动学方程+欧拉参数方程+移动位移方程)

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

if nargin < 3 || isempty(flag)
    
    %% 由积分变量增加值更新各矩阵（广义速率、欧拉参数、移动位移）
    ntemp=size(Var,2);  
    % 更新未知广义速率阵列
    y(Var)=col(1:ntemp);  
    % 更新欧拉参数矩阵
    epsi=reshape(col(ntemp+1:ntemp+N*4),4,N);
    % 更新相对位移矢量矩阵
    for i=1:length(Varmove_row)
        s(Varmove_row(i),Varmove_col(i))=col(ntemp+4*N+i);
    end

    %% ★计算相对变换矩阵及绝对变换矩阵★(相对变换矩阵与欧拉参数有关)
    % 注意此时的计算方式与主程序中计算SJK0不同
    % 计算相对变换矩阵[SJK]
    % SJK(:,:,k)为用欧拉参数表示的刚体Bk向其内接刚体Bj转换的变换矩阵
    E=eye(3);
    for k=1:N
        for i=1:3
            for j=1:3
                SJKtemp(i,j)=epsi(i,k)*epsi(j,k);
            end
        end
        SJK(:,:,k)=2*(SJKtemp+epsi(4,k)*(epsi(4,k)*E+Dual_Mat(epsi(1:3,k))))-E;
    end
    % 计算绝对变换矩阵[SOK]
    % SOK(:,:,k)为用欧拉参数表示的刚体Bk向静止坐标系转换的变换矩阵
    for k=1:N
        SOK(1:3,1:3,k)=E;
        for i=PathN(k):-1:1
            SOK(:,:,k)=SOK(:,:,k)*SJK(:,:,Path(i,k));
        end
    end
    
    %% ★计算偏角速度、偏速度及它们的导数★（与广义速率、绝对变换矩阵、绝对变换矩阵的导数、连接点矢量、质心矢量、移动变量、移动速率有关）
    % Pomiga偏角速度在静止坐标系中的投影分量矩阵(即公式中的Omiga_kml)
    Pomiga(1:Free,1:3,1:N)=0;
    Pomiga(1:3,:,1)=E;
    for k=2:N
        Pomiga(1:3*k-3,:,k)=Pomiga(1:3*k-3,:,Low(k));
        % SOK.'是因为Omiga_kml要与Smp相对应
        Pomiga(3*k-2:3*k,:,k)=SOK(:,:,Low(k)).';
    end
    % 绝对变换矩阵导数SOKdot
    for k=1:N
        % Omiga(:,k)表示刚体Bk的角速度在静止坐标系中的投影分量矩阵(按列向量记录)
        Omiga(:,k)=Pomiga(1:Free,:,k)'*y(1:Free);
        SOKdot(:,:,k)=Dual_Mat(Omiga(:,k))*SOK(:,:,k);
    end
    % 偏角速度的导数Pomigadot
    Pomigadot(1:Free,1:3,1:N)=0;
    for k=2:N
        Pomigadot(1:3*k-3,:,k)=Pomigadot(1:3*k-3,:,Low(k));
        Pomigadot(3*k-2:3*k,:,k)=SOKdot(:,:,Low(k)).';
    end
    % 偏速度Pv
    Pv(1:Free,1:3,1:N)=0;
    for k=1:N
        for i=1:PathN(k)-1
            I=Path(i+1,k);  J=Path(i,k);     % 用I来表示J的内接刚体
            Pv(1:3*k,:,k)=Pv(1:3*k,:,k)...
                +Pomiga(1:3*k,:,I)*Dual_Mat(SOK(:,:,I)*(Q_Local(:,J)+s(:,J)));
        end
        Pv(1:3*k,:,k)=Pv(1:3*k,:,k)...
            +Pomiga(1:3*k,:,k)*Dual_Mat(SOK(:,:,k)*R_Local(:,k));
        Pv(3*N+1:6*N,:,1:N)=Pomiga(1:3*N,:,1:N);
    end
    % 偏速度的导数Pvdot
    Pvdot(1:Free,1:3,1:N)=0;
    for k=1:N
        for i=1:PathN(k)-1
            I=Path(i+1,k);  J=Path(i,k);
            sdot(:,J)=y(3*(N+J-1)+1:3*(N+J-1)+3);
            Pvdot(1:3*k,:,k)=Pvdot(1:3*k,:,k)...
                +Pomigadot(1:3*k,:,I)*Dual_Mat(SOK(:,:,I)*(Q_Local(:,J)+s(:,J)))...
                +Pomiga(1:3*k,:,I)*Dual_Mat(SOKdot(:,:,I)*(Q_Local(:,J)+s(:,J)))...
                +Pomiga(1:3*k,:,I)*Dual_Mat(SOK(:,:,I)*sdot(:,J));
        end
        Pvdot(1:3*k,:,k)=Pvdot(1:3*k,:,k)...
            +Pomigadot(1:3*k,:,k)*Dual_Mat(SOK(:,:,k)*R_Local(:,k))...
            +Pomiga(1:3*k,:,k)*Dual_Mat(SOKdot(:,:,k)*R_Local(:,k));
        Pvdot(3*N+1:6*N,:,1:N)=Pomigadot(1:3*N,:,1:N);
    end
    
    %% ★生成基本运动力学方程系数阵a_kl和f_1★
    %% 计算广义主动力和广义惯性力等(与Dyn_Res_of_Anc_Cha_Anc1_Fun_Bk中计算重复)
    % 广义主动力F_Act<6N*1>为列向量
    F_Act(1:Free,1)=0;  
    nt=size(Br,2);
    F_Act=F_Act+Pv(:,:,N)*FR(1:3,nt)+Pomiga(:,:,N)*FR(4:6,nt)...
        +Pv(:,:,N)*FC(1:3,nt,N);

    F_Act=F_Act+Pv(:,:,1)*FF(1:3,nt);  
  
    for k=1:N
        F_Act=F_Act+Pv(:,1,k)*(F_Vx(nt,k)+F_Ax(nt,k)+F_Dx(nt,k))...
            +Pv(:,3,k)*(F_Vz(nt,k)+F_Az(nt,k)+F_Dz(nt,k))...
            +Pomiga(:,2,k)*(T_A(nt,k)+T_D(nt,k)+T_V(nt,k))...
            +Pv(:,:,k)*(F(:,k)+FB(:,k)+FSC(1:3,nt,k)+SOK(:,:,k)*F_Local(:,k))...
            +Pomiga(:,:,k)*(M(:,k)+MSC(1:3,nt,k)+SOK(:,:,k)*M_Local(:,k)); %包括局部力和整体力
    end

    % 通过坐标变换计算静止坐标系中的中心惯性矩
    for k=1:N
        % 为什么是这样的变换公式???《浮-沈》p26
        IT(:,:,k)=SOK(:,:,k)'*IT_Local(:,:,k)*SOK(:,:,k);
    end
    % 计算广义惯性力
    % a_lp为广义惯性力中与广义速率导数(加速度)有关的系数矩阵,即a_lp
    % h_l为广义惯性力中与广义速率有关的项,即h_l
    a_lp(1:Free,1:Free)=0;  h_l(1:Free,1)=0;     %对称且正定
    for k=1:N
        a_lp=a_lp+Mass(k)*Pv(:,:,k)*Pv(:,:,k)'...
            +Pomiga(:,:,k)*IT(:,:,k)*Pomiga(:,:,k)';
        h_l=h_l+Mass(k)*Pv(:,:,k)*Pvdot(:,:,k)'*y...
            +Pomiga(:,:,k)*IT(:,:,k)*Pomigadot(:,:,k)'*y...
            +Pomiga(:,:,k)*Dual_Mat(Pomiga(:,:,k)'*y)*IT(:,:,k)*(Pomiga(:,:,k)'*y);
    end
    

    %% ★构建约束方程★
    % 二维锚链多体系统,其约束方程可写为B_Sub*y(Var)=g_Ext或Bsub*Ysub=g_Ext;Ysub=y(Var)
    % g_Ext_xz表示x、z方向的外部约束，表示z方向的外部约束
    g_Ext_xz=zeros(2,1);      g_Ext_z=zeros(1,1);
    g_Extdot_xz=zeros(2,1);   g_Extdot_z=zeros(1,1);
    % 求解外部约束矩阵B_Ext及其导数Bdot
    B_Ext1=Pomiga(:,:,1).';  B_Ext2=Pv(:,:,1).';  B_Ext=[B_Ext1; B_Ext2];
    % 求解缩减后的约束矩阵Bsub[m*r];  r=length(Var)
    Bsub_xz=B_Ext([4 6],Var);   Bsub_z=B_Ext(6,Var);
    % 求解Bsub的正交补阵Csub[r*(r-m)];
    % Method 1 - Csub is a matrix composed of r-m eigenvectors which ...
    % correspond to the zero eigenvalues of Bsub.'*Bsub.
    % 方法1-Csub是以(r-m)个Bsub.'*Bsub的零特征跟相应的特征矢量为列的矩阵。
    [m_xz,r_xz]=size(Bsub_xz);  [C1_xz,D_xz]=eig(Bsub_xz.'*Bsub_xz);      Csub_xz=[];
    % 选出Bsub.'*Bsub零特征根对应的特征向量，组成矩阵C
    for k=1:r_xz
        if abs(D_xz(k,k))<=1.0e-6
            Csub_xz=[Csub_xz  C1_xz(:,k)];
        end
    end
    [m_z,r_z]=size(Bsub_z);  [C1_z,D_z]=eig(Bsub_z.'*Bsub_z);      Csub_z=[];
    for k=1:r_z
        if abs(D_z(k,k))<=1.0e-6
            Csub_z=[Csub_z  C1_z(:,k)];
        end
    end
    % 求解Bsub的导数Bsubdot
    B_Ext1dot=Pomigadot(:,:,1).';     B_Ext2dot=Pvdot(:,:,1).';
    B_Extdot=[B_Ext1dot; B_Ext2dot];
    Bsubdot_xz=B_Extdot([4 6],Var);
    Bsubdot_z=B_Extdot(6,Var);
    % 系统约束方程,两边求导后为Bsub*ydot(Var)=g_Extdot-Bsubdot*y(Var)
    
    %% 组集运动方程，缩减自由度
    % Fsum为与广义速率有关的广义力
    % Asub为缩减自由度后的系统基本运动微分方程的左端系数矩阵
    % inv(Asub)*Fsub为待求解的微分方程的右边部分
    % 缩减后的运动力学方程为Asub*ydot(Var)=Fsub+Bsub.'*lamda...
    % 两边同乘以Csub.'可得Csub.'*Asub*ydot(Var)=Csub.'*Fsub
    Fsum=F_Act-h_l;
    Asub=a_lp(Var,Var);
    Fsub=Fsum(Var);
    % 系统约束方程与缩减后的系统运动力学方程联立后变为UA*ydot(Var)=UB
    UA_xz=[Csub_xz.'*Asub;  Bsub_xz];
    UB_xz=[Csub_xz.'*Fsub;  g_Extdot_xz-Bsubdot_xz*y(Var)];
    UA_z=[Csub_z.'*Asub;  Bsub_z];
    UB_z=[Csub_z.'*Fsub;  g_Extdot_z-Bsubdot_z*y(Var)];
    % 作为下一时刻待求自由度的广义速率的导数
    % F_Hin_Kane为负时表示拉力，正时表示压力
    if (nt>=2&&-F_Hin_Kane(1,nt-1,2)>-FDF(1,nt))||Bv(1,nt,1)>1.0e-6
        % 作为下一时刻待求自由度的广义速率的导数
        dy=UA_z\UB_z;
    else
        % 该方程需要如何处理？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
        dy=UA_xz\UB_xz;  
    end
    
    % depsi为欧拉参数变换矩阵的导数(与欧拉参数和广义速率有关)
    for k=1:N
        depsi(4,k)=-epsi(1:3,k)'*y(3*k-2:3*k)/2;
        depsi(1:3,k)=(-Dual_Mat(epsi(1:3,k))+epsi(4,k)*E)*y(3*k-2:3*k)/2;
    end
    ds=y(Varmove);
    % 定义输出的参数,若不缩减则输出参数的个数为(6N+4N+3N)
    out1=[dy; reshape(depsi,4*N,1); ds];
    
else
    
    switch(flag)
        case 'init'
            % 初始化设置,此时调用[tspan,y0,options]=odefun([],[],'init')
            out1=0:Delta_t:Step_n*Delta_t;
            out2=col0;
            out3=odeset('reltol',0.001,'maxstep',Delta_t,'InitialStep',Delta_t/2);
        otherwise
            error(['Unknown flag ''' flag '''.']);
    end
end