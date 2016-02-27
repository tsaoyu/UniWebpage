function Dyn_Res_of_Anc_Cha_Anc1_Fun_Bk(col,nt)
% nt=1～(Step_n+1)为积分时刻序号,
% col为一列向量,其元素为[待求自由度广义速率;欧拉参数;移动位移]
% 本函数的功能为输入nt时刻及此时刻微分方程方程的解,从而通过Kane方法计算...
% 系统的动力响应(不同时刻刚体Bk的物理量及广义速率的导数+接点约束力)

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
global dK_nt  E_Flyl  E_Kin_k  E_Gra_k  E_Kin  E_Gra  E0_Kin
global mu  FR  E_FR  FDF  FF  FSC
% 基于牛顿定律、Kane方法求解接点约束力的相关变量
% T_Hin为通过Kane方法计算得到的物体所受力矩：M_After为经典力学方法计算得到的物体所受力矩
% M_Imp是使锚链平衡所施加的力矩
global T_Hin  M_After  MSC 
% F_Con为广义约束力,F_Hin_Kane为Kane方法计算得到的接点约束力
global F_Con  F_Hin_Kane  Mcon_3N  Fcon_3N  Mcon_k

% ydottemp在主程序中没有定义,其主要功能是将nt时刻待求自由度所对应广义速率的
% 导数ydot传递给(nt+1)时刻.这样做可以实现！
% 因为在本函数中声明的全局变量可以保存在本函数的工作空间中,若不声明则变量会被清除！
% 当然,在主程序中加上定义global yottemp也可以.
global ydottemp 

%　附加质量力附加阻尼力及粘性摩擦阻力相关
global  v_t  v_C  v_Cdot  F_V  T_V  F_A  T_A  F_D  T_D   F_Vx   F_Vz  F_Ax   F_Az  F_Dx   F_Dz 
% 外部约束相关变量
global  gt_xz  gtdot_xz  Error_e_xz  lamda_xz  gt_z  gtdot_z  Error_e_z  lamda_z  Bsub_xz  Csub_xz  Bsub_z  Csub_z 

%% 读取广义速率的导数 
% 初始时刻nt=1，Br=[]求Br;ydot(Var)已知.之后Br不为空
if ~isempty(Br)
   if nt>size(Br,2) 		%在每一 nt 时刻只计算一次,此时nt为新一时刻 
      ydot(Var)=ydottemp;
   else                     %此时nt=size(Br,2)说明已经计算过
      return;
   end
end
    
if nt/100==fix(nt/100)
    disp(['nt =' num2str(nt) ]); 
end


%% 更新各矩阵（广义速率、欧拉参数、移动位移）
% Var为待求解自由度的编号(是一行向量)
% ntemp为待求解自由度的数量
ntemp=size(Var,2);
% 更新未知广义速率阵列
y(Var)=col(1:ntemp);
% 更新欧拉参数矩阵
epsi=reshape(col(ntemp+1:ntemp+N*4),4,N);
% 更新相对位移矢量矩阵
for i=1:length(Varmove_row)
    s(Varmove_row(i),Varmove_col(i))=col(ntemp+4*N+i);
end

%% ★计算相对变换矩阵及绝对变换矩阵★
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

%% ★计算偏角速度、偏速度及它们的导数★
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
    Omiga(:,k)=Pomiga(1:3*k,:,k)'*y(1:3*k);
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
        sdot(:,J)=y(3*N+3*(J-1)+1:3*N+3*(J-1)+3);
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

%% 求各时刻刚体的原点坐标和质心坐标
% p0(:,k)为刚体Bk的连体坐标系原点在静止坐标系中的投影
p0(:,1)=Q_Local(:,1)+s(:,1);
for k=2:N
    p0(:,k)=p0(:,Low(k))+SOK(:,:,Low(k))*(Q_Local(:,k)+s(:,k));
end
% pG(:,k)为刚体Bk质心Gk在静止坐标系中的坐标阵列
% PV0(j,i,k)为刚体Bk八个顶点相对其质心的xyz坐标(连体坐标系中)
% PV(:,i,k,nt)为刚体Bk的八个顶点在nt时刻的xyz坐标(在静止坐标系中)
for k=1:N
    pG(:,k)=p0(:,k)+SOK(:,:,k)*R_Local(:,k);
    for i=1:8
        PV(:,i,k,nt)=pG(:,k)+(SOK(:,:,k)*PV0(:,i,k));
    end
end

%% ★在低序体连体坐标系内求相对运动响应★
% Brj(1:3,nt,k)表示刚体Bk在nt时刻相对于Bj的角位移,其他变量含义类似,一般无需计算
% 注意角度(钝角)的判断
for k=1:N
    Brj(2,nt,k)=real(asin(SJK(1,3,k)));
    Brj(1,nt,k)=real(asin(-SJK(2,3,k)/cos(Brj(2,nt,k))));
    Brj(3,nt,k)=real(asin(-SJK(1,2,k)/cos(Brj(2,nt,k))));
    Bwj(:,nt,k)=y(3*k-2:3*k);
    Bej(:,nt,k)=ydot(3*k-2:3*k);    
    
    Bsj(:,nt,k)=s(:,k);
    Bvj(:,nt,k)=y(3*(k+N-1)+1:3*(k+N-1)+3);
    Baj(:,nt,k)=ydot(3*(k+N-1)+1:3*(k+N-1)+3);
end

%% ★在惯性坐标系内求绝对运动响应★
% Br(1:3,1:stepN,1:N)表示各刚体相对静止坐标系的角位移
% 注意角度(-pi～pi)的判断,
% Bw(1:3,1:stepN,1:N)表示各刚体的角速度
% Be(1:3,1:stepN,1:N)表示各刚体的角加速度
% Bs(1:3,1:stepN,1:N)表示各刚体质心的线位移
% Bv(1:3,1:stepN,1:N)表示各刚体质心的线速度
% Ba(1:3,1:stepN,1:N)表示各刚体质心的线加速度
for k=1:N

    Br(2,nt,k)=real(asin(SOK(1,3,k)));
    Br(1,nt,k)=real(asin(-SOK(2,3,k)/cos(Br(2,nt,k))));
    Br(3,nt,k)=real(asin(-SOK(1,2,k)/cos(Br(2,nt,k))));
    Bs(:,nt,k)=pG(:,k);
        
    % Angle为与x0y0z0轴的正向夹角
    Angle(:,nt,k)=-Br(:,nt,k);
    
    Bw(:,nt,k)=Pomiga(:,:,k)'*y;
    Bv(:,nt,k)=Pv(:,:,k)'*y;
    Be(:,nt,k)=Pomigadot(:,:,k)'*y+Pomiga(:,:,k)'*ydot;
    Ba(:,nt,k)=Pvdot(:,:,k)'*y+Pv(:,:,k)'*ydot;
      
    % 通过坐标变换计算静止坐标系中的中心惯性矩
    IT(:,:,k)=SOK(:,:,k)'*IT_Local(:,:,k)*SOK(:,:,k);
 
    % 计算nt时刻,系统的动能(质心的平动+相对质心平移坐标系的转动)
    E_Kin_k(k,nt)=Mass(k)*sum(Bv(:,nt,k).^2)/2+...
        diag(IT_Local(:,:,k)).'*Bw(:,nt,k).^2/2;
     
    % 计算nt时刻,系统的势能(重力-浮力)★★★注意不同
    E_Gra_k(k,nt)=(-Mass(k)*g-FB(3,k))*(pG(3,k)+h_Pre+L_Aux+R_Cha);  
    
end
E_Kin(nt)=sum(E_Kin_k(:,nt));
E_Gra(nt)=sum(E_Gra_k(:,nt));
% 船舶停止标志位
if Bv(1,nt,N)<=1.0e-6
    BvSFlag=1;
end

%% ★生成基本运动力学方程系数阵a_kl和f_1★
%% 计算广义主动力,广义惯性力,并进行数值结果验证(计算nt时刻的受力供nt+1时刻用及ydottemp)
% 广义主动力F_Act<6N*1>为列向量
F_Act(1:Free,1)=0;  
nt=size(Br,2);
% 因为小角度时sin(theta)=theta???
% =========================================================================  ===================== 待改进 ==========================
base_dis=[Bs(:,nt,N);sin(Br(:,nt,N))];
base_acc=[Ba(:,:,N)',Be(:,:,N)']; 
% =========================================================================  ===================== 待改进 ==========================
% FC为水流力矩阵
vS(1,nt)=Bv(1,nt,N);           vW(1,nt)=1.261;                                                 
vSW(1,nt)=vS(1,nt)-vW(1,nt);    
Re=abs(vSW(1,nt))*LenN(1)/upsilon;
Cf=0.075/(log(Re)-2)^2;
Csf=(105*(ks/LenN(1))^(1/3)-0.64)*1.0e-3; 
FC(1,nt,N)=-sign(vSW(1,nt))*(Cf+Csf)*rou_W*S_Wet*vSW(1,nt).^2/2; 
FC(2,nt,N)=-sign(vSW(2,nt))*rou_W*Cy*Ay*vSW(2,nt).^2/2; 
% FR为静水恢复力矩阵；
FR(:,nt)=FR_Factor.*base_dis;
E_FR(nt)=-sum(FR_Factor.*[Bs(:,nt,N).^2;Br(:,nt,N).^2])/2;
F_Act=F_Act+Pv(:,:,N)*FR(1:3,nt)+Pomiga(:,:,N)*FR(4:6,nt)...
    +Pv(:,:,N)*FC(1:3,nt,N);

% =========================================================================  ==================== 待改进 ==========================
% 取消弹性地基梁对重力锚的作用，而用外力约束代替
for k=2:N-2
    if  Bv(3,nt,k)<0&&Bs(3,nt,k)<-h_Pre-L_Aux
        if k<=N_Cha+1
            FSC(3,nt,k)=-c_Seabed*(Bs(3,nt,k)+h_Pre+L_Aux)*prod(Len_Cha(1:2));
        else
            FSC(3,nt,k)=-c_Seabed*(Bs(3,nt,k)+h_Pre+L_Aux)*prod(Len_Pre(1:2));
        end
    end
end

% =========================================================================  ==================== 待改进 ==========================
FDF(1,1)=-mu*(-Mass(1)*g-FB(3,1));  
if nt>=2
    % F_Hin_Kane为负值时表示拉力
    FDF(1:3,nt)=[-mu*(-Mass(1)*g-FB(3,1)+F_Hin_Kane(3,nt-1,2)); 0; 0];
end
if Bv(1,nt,1)>0
    FF(1,nt)=FDF(1,nt); 
end
F_Act=F_Act+Pv(:,:,1)*FF(1:3,nt);

for k=2:N-2
    v_t(nt,k)=-Bv(1,nt,k)*cos(-Br(2,nt,k)) - Bv(3,nt,k)*sin(-Br(2,nt,k));    
    v_Cdot(nt,k)=Ba(1,nt,k)*sin(-Br(2,nt,k)) - Ba(3,nt,k)*cos(-Br(2,nt,k));   
    v_C(nt,k)=Bv(1,nt,k)*sin(-Br(2,nt,k)) - Bv(3,nt,k)*cos(-Br(2,nt,k));    
    F_V(nt,k)=sign(v_t(nt,k))*pi*rou_W*(2*Rad(k))*Len(1,k)*C_F*v_t(nt,k)^2 / 2;
    T_V(nt,k)=0;
    F_A(nt,k)=pi*rou_W*(2*Rad(k))^2*C_M*v_Cdot(nt,k)*Len(1,k) / 4;
    T_A(nt,k)=-pi*rou_W*(2*Rad(k))^2*C_M*Be(2,nt,k)*Len(1,k)^3 / 48;
    if abs(v_C(nt,k))>=abs(Bw(2,nt,k))*Len(1,k)/2
        F_D(nt,k)=sign(v_C(nt,k))*rou_W*(2*Rad(k))*C_D* ( v_C(nt,k)^2*Len(1,k) + Bw(2,nt,k)^2*Len(1,k)^3/12 ) / 2;
        T_D(nt,k)=-sign(Bw(2,nt,k))*rou_W*(2*Rad(k))*C_D*abs(v_C(nt,k)*Bw(2,nt,k))*Len(1,k)^3 / 12;
    elseif abs(v_C(nt,k))<abs(Bw(2,nt,k))*Len(1,k)/2
        F_D(nt,k)=rou_W*(2*Rad(k))*C_D* ( 2*v_C(nt,k)^3/( 3*abs(Bw(2,nt,k))) + v_C(nt,k)*abs(Bw(2,nt,k))*Len(1,k)^2/2 ) / 2;
        T_D(nt,k)=-sign(Bw(2,nt,k))*rou_W*(2*Rad(k))*C_D* ( Bw(2,nt,k)^2*Len(1,k)^4/32 + v_C(nt,k)^2*Len(1,k)^2/4 - v_C(nt,k)^4/(6*Bw(2,nt,k)^2) ) / 2 ;
    end
        
    F_Vx(nt,k)=F_V(nt,k)*sin(-Br(2,nt,k));   F_Vz(nt,k)=F_V(nt,k)*cos(-Br(2,nt,k));
    F_Ax(nt,k)=-F_A(nt,k)*sin(-Br(2,nt,k));  F_Az(nt,k)=F_A(nt,k)*cos(-Br(2,nt,k));
    F_Dx(nt,k)=-F_D(nt,k)*sin(-Br(2,nt,k));  F_Dz(nt,k)=F_D(nt,k)*cos(-Br(2,nt,k));
       
end 
 
for k=1:N
    F_Act=F_Act+Pv(:,1,k)*(F_Vx(nt,k)+F_Ax(nt,k)+F_Dx(nt,k))...
               +Pv(:,3,k)*(F_Vz(nt,k)+F_Az(nt,k)+F_Dz(nt,k))...
               +Pomiga(:,2,k)*(T_A(nt,k)+T_D(nt,k)+T_V(nt,k))...
               +Pv(:,:,k)*(F(:,k)+FB(:,k)+FSC(:,nt,k)+SOK(:,:,k)*F_Local(:,k))...
               +Pomiga(:,:,k)*(M(:,k)+MSC(1:3,nt,k)+SOK(:,:,k)*M_Local(:,k)); %包括局部力和整体力
end

% nt时刻,广义主动力F_Act与广义速率y的乘积(对应元素相乘后求和).nt=1:Step_n+1
% 对于某些非完整系统，需要加常数项！☆☆☆
dK_nt(nt)=Delta_t*F_Act.'*y;
E_Flyl(nt)=sum(dK_nt(1:nt))+E0_Kin;
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
    ydottemp=UA_z\UB_z;
    gt_z(:,nt)=Bsub_z*y(Var);
    gtdot_z(:,nt)=Bsub_z*ydottemp+Bsubdot_z*y(Var);
    Error_e1_z(:,nt)=Csub_z.'*Asub*ydottemp-Csub_z.'*Fsub;
    Error_e2_z(:,nt)=Bsub_z*ydottemp-g_Extdot_z+Bsubdot_z*y(Var);
    Error_e_z=[Error_e1_z;  Error_e2_z];
    % 求解约束力和约束力矩的诸分量lamda
    lamda_z(:,nt)=Bsub_xz.'\(Asub*ydottemp-Fsub);
else
    % 该方程需要如何处理？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
    ydottemp=UA_xz\UB_xz;
    % gt为记录不同时刻的约束； yz为方程Bsub*ydottemp-g_Extdot+Bsubdot*y(Var)的误差
    gt_xz(:,nt)=Bsub_xz*y(Var);
    gtdot_xz(:,nt)=Bsub_xz*ydottemp+Bsubdot_xz*y(Var);
    Error_e1_xz(:,nt)=Csub_xz.'*Asub*ydottemp-Csub_xz.'*Fsub;
    Error_e2_xz(:,nt)=Bsub_xz*ydottemp-g_Extdot_xz+Bsubdot_xz*y(Var);
    Error_e_xz=[Error_e1_xz;  Error_e2_xz];
    % 求解约束力和约束力矩的诸分量lamda
    lamda_xz(:,nt)=Bsub_xz.'\(Asub*ydottemp-Fsub);
end

%% 基于Kane方法求解接点约束力
% F_con为广义约束力矩阵，经过相关变换后可求得各接点处的约束力
Arem=a_lp(Var_Know,Var);
Frem=Fsum(Var_Know);
% F_Con为广义约束力！！！
F_Con(1:Free,nt)=0; 
F_Con(Var_Know,nt)=Arem*ydottemp-Frem;
% Mcon_3N和Fcon_3N分别为相对于各物体质心处所施加的约束力矩和约束力（非广义力）
% 经过变换，可求得其连接点处的约束力力矩和约束力
% Pomiga_Reshape*Mcon_3N+Pv_Reshape*Fcon_3N=F_Con(大小为6N*1); 
Pomiga_Reshape(:,:)=reshape(Pomiga,6*N,3*N);
Pv_Reshape(:,:)=reshape(Pv,6*N,3*N);  
P_6N=[Pomiga_Reshape Pv_Reshape];
Fcon_6N(:,nt)=P_6N\F_Con(:,nt);
Mcon_3N(:,nt)=Fcon_6N(1:3*N,nt);  Fcon_3N(:,nt)=Fcon_6N(3*N+1:6*N,nt);
% 此多级摆系统中，Mcon_k等于M_Hin
Mcon_k(:,nt,:)=reshape(Mcon_3N(:,nt),3,N);
Fcon_k(:,nt,:)=reshape(Fcon_3N(:,nt),3,N);
for k=1:N
    F_Hin_Kane(:,nt,k)=sum(Fcon_k(:,nt,k:N),3); 
end

% 验证海底约束力与约束力矩
for k=2:N-2
    M_After(:,nt,k)=diag(IT_Local(:,:,k)).*Be(:,nt,k);
    if Angle(2,nt,k)>=-pi/2&&Angle(2,nt,k)<0
        T_Hin(2,nt,k)=(F_Hin_Kane(1,nt,k)+F_Hin_Kane(1,nt,k+1))*Len(1,k)/2*sin(-Angle(2,nt,k))+...
            (F_Hin_Kane(3,nt,k)+F_Hin_Kane(3,nt,k+1))*Len(1,k)/2*cos(-Angle(2,nt,k));
    elseif Angle(2,nt,k)>=0&&Angle(2,nt,k)<=pi/2
        T_Hin(2,nt,k)=(-F_Hin_Kane(1,nt,k)-F_Hin_Kane(1,nt,k+1))*Len(1,k)/2*sin(Angle(2,nt,k))+...
            (F_Hin_Kane(3,nt,k)+F_Hin_Kane(3,nt,k+1))*Len(1,k)/2*cos(Angle(2,nt,k));
    end

end

