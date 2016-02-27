function Dyn_Res_of_Anc_Cha_Anc1_Fun_Bk(col,nt)
% nt=1��(Step_n+1)Ϊ����ʱ�����,
% colΪһ������,��Ԫ��Ϊ[�������ɶȹ�������;ŷ������;�ƶ�λ��]
% �������Ĺ���Ϊ����ntʱ�̼���ʱ��΢�ַ��̷��̵Ľ�,�Ӷ�ͨ��Kane��������...
% ϵͳ�Ķ�����Ӧ(��ͬʱ�̸���Bk�����������������ʵĵ���+�ӵ�Լ����)

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
global dK_nt  E_Flyl  E_Kin_k  E_Gra_k  E_Kin  E_Gra  E0_Kin
global mu  FR  E_FR  FDF  FF  FSC
% ����ţ�ٶ��ɡ�Kane�������ӵ�Լ��������ر���
% T_HinΪͨ��Kane��������õ��������������أ�M_AfterΪ������ѧ��������õ���������������
% M_Imp��ʹê��ƽ����ʩ�ӵ�����
global T_Hin  M_After  MSC 
% F_ConΪ����Լ����,F_Hin_KaneΪKane��������õ��Ľӵ�Լ����
global F_Con  F_Hin_Kane  Mcon_3N  Fcon_3N  Mcon_k

% ydottemp����������û�ж���,����Ҫ�����ǽ�ntʱ�̴������ɶ�����Ӧ�������ʵ�
% ����ydot���ݸ�(nt+1)ʱ��.����������ʵ�֣�
% ��Ϊ�ڱ�������������ȫ�ֱ������Ա����ڱ������Ĺ����ռ���,��������������ᱻ�����
% ��Ȼ,���������м��϶���global yottempҲ����.
global ydottemp 

%������������������������ճ��Ħ���������
global  v_t  v_C  v_Cdot  F_V  T_V  F_A  T_A  F_D  T_D   F_Vx   F_Vz  F_Ax   F_Az  F_Dx   F_Dz 
% �ⲿԼ����ر���
global  gt_xz  gtdot_xz  Error_e_xz  lamda_xz  gt_z  gtdot_z  Error_e_z  lamda_z  Bsub_xz  Csub_xz  Bsub_z  Csub_z 

%% ��ȡ�������ʵĵ��� 
% ��ʼʱ��nt=1��Br=[]��Br;ydot(Var)��֪.֮��Br��Ϊ��
if ~isempty(Br)
   if nt>size(Br,2) 		%��ÿһ nt ʱ��ֻ����һ��,��ʱntΪ��һʱ�� 
      ydot(Var)=ydottemp;
   else                     %��ʱnt=size(Br,2)˵���Ѿ������
      return;
   end
end
    
if nt/100==fix(nt/100)
    disp(['nt =' num2str(nt) ]); 
end


%% ���¸����󣨹������ʡ�ŷ���������ƶ�λ�ƣ�
% VarΪ��������ɶȵı��(��һ������)
% ntempΪ��������ɶȵ�����
ntemp=size(Var,2);
% ����δ֪������������
y(Var)=col(1:ntemp);
% ����ŷ����������
epsi=reshape(col(ntemp+1:ntemp+N*4),4,N);
% �������λ��ʸ������
for i=1:length(Varmove_row)
    s(Varmove_row(i),Varmove_col(i))=col(ntemp+4*N+i);
end

%% �������Ա任���󼰾��Ա任�����
% ������Ա任����[SJK]
% SJK(:,:,k)Ϊ��ŷ��������ʾ�ĸ���Bk�����ڽӸ���Bjת���ı任����
E=eye(3);
for k=1:N
    for i=1:3
        for j=1:3
            SJKtemp(i,j)=epsi(i,k)*epsi(j,k);
        end
    end
    SJK(:,:,k)=2*(SJKtemp+epsi(4,k)*(epsi(4,k)*E+Dual_Mat(epsi(1:3,k))))-E;
end
% ������Ա任����[SOK]
% SOK(:,:,k)Ϊ��ŷ��������ʾ�ĸ���Bk��ֹ����ϵת���ı任����
for k=1:N
    SOK(1:3,1:3,k)=E;
    for i=PathN(k):-1:1
        SOK(:,:,k)=SOK(:,:,k)*SJK(:,:,Path(i,k));
    end
end

%% �����ƫ���ٶȡ�ƫ�ٶȼ����ǵĵ�����
% Pomigaƫ���ٶ��ھ�ֹ����ϵ�е�ͶӰ��������(����ʽ�е�Omiga_kml)
Pomiga(1:Free,1:3,1:N)=0;
Pomiga(1:3,:,1)=E;
for k=2:N
    Pomiga(1:3*k-3,:,k)=Pomiga(1:3*k-3,:,Low(k));
    % SOK.'����ΪOmiga_kmlҪ��Smp���Ӧ
    Pomiga(3*k-2:3*k,:,k)=SOK(:,:,Low(k)).';
end
% ���Ա任������SOKdot
for k=1:N
    % Omiga(:,k)��ʾ����Bk�Ľ��ٶ��ھ�ֹ����ϵ�е�ͶӰ��������(����������¼)
    Omiga(:,k)=Pomiga(1:3*k,:,k)'*y(1:3*k);
    SOKdot(:,:,k)=Dual_Mat(Omiga(:,k))*SOK(:,:,k);
end
% ƫ���ٶȵĵ���Pomigadot
Pomigadot(1:Free,1:3,1:N)=0;
for k=2:N
    Pomigadot(1:3*k-3,:,k)=Pomigadot(1:3*k-3,:,Low(k));
    Pomigadot(3*k-2:3*k,:,k)=SOKdot(:,:,Low(k)).';
end
% ƫ�ٶ�Pv
Pv(1:Free,1:3,1:N)=0;
for k=1:N
    for i=1:PathN(k)-1
        I=Path(i+1,k);  J=Path(i,k);     % ��I����ʾJ���ڽӸ���
        Pv(1:3*k,:,k)=Pv(1:3*k,:,k)...
            +Pomiga(1:3*k,:,I)*Dual_Mat(SOK(:,:,I)*(Q_Local(:,J)+s(:,J)));
    end
    Pv(1:3*k,:,k)=Pv(1:3*k,:,k)...
        +Pomiga(1:3*k,:,k)*Dual_Mat(SOK(:,:,k)*R_Local(:,k));
    Pv(3*N+1:6*N,:,1:N)=Pomiga(1:3*N,:,1:N);
end
% ƫ�ٶȵĵ���Pvdot
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

%% ���ʱ�̸����ԭ���������������
% p0(:,k)Ϊ����Bk����������ϵԭ���ھ�ֹ����ϵ�е�ͶӰ
p0(:,1)=Q_Local(:,1)+s(:,1);
for k=2:N
    p0(:,k)=p0(:,Low(k))+SOK(:,:,Low(k))*(Q_Local(:,k)+s(:,k));
end
% pG(:,k)Ϊ����Bk����Gk�ھ�ֹ����ϵ�е���������
% PV0(j,i,k)Ϊ����Bk�˸�������������ĵ�xyz����(��������ϵ��)
% PV(:,i,k,nt)Ϊ����Bk�İ˸�������ntʱ�̵�xyz����(�ھ�ֹ����ϵ��)
for k=1:N
    pG(:,k)=p0(:,k)+SOK(:,:,k)*R_Local(:,k);
    for i=1:8
        PV(:,i,k,nt)=pG(:,k)+(SOK(:,:,k)*PV0(:,i,k));
    end
end

%% ���ڵ�������������ϵ��������˶���Ӧ��
% Brj(1:3,nt,k)��ʾ����Bk��ntʱ�������Bj�Ľ�λ��,����������������,һ���������
% ע��Ƕ�(�۽�)���ж�
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

%% ���ڹ�������ϵ��������˶���Ӧ��
% Br(1:3,1:stepN,1:N)��ʾ��������Ծ�ֹ����ϵ�Ľ�λ��
% ע��Ƕ�(-pi��pi)���ж�,
% Bw(1:3,1:stepN,1:N)��ʾ������Ľ��ٶ�
% Be(1:3,1:stepN,1:N)��ʾ������ĽǼ��ٶ�
% Bs(1:3,1:stepN,1:N)��ʾ���������ĵ���λ��
% Bv(1:3,1:stepN,1:N)��ʾ���������ĵ����ٶ�
% Ba(1:3,1:stepN,1:N)��ʾ���������ĵ��߼��ٶ�
for k=1:N

    Br(2,nt,k)=real(asin(SOK(1,3,k)));
    Br(1,nt,k)=real(asin(-SOK(2,3,k)/cos(Br(2,nt,k))));
    Br(3,nt,k)=real(asin(-SOK(1,2,k)/cos(Br(2,nt,k))));
    Bs(:,nt,k)=pG(:,k);
        
    % AngleΪ��x0y0z0�������н�
    Angle(:,nt,k)=-Br(:,nt,k);
    
    Bw(:,nt,k)=Pomiga(:,:,k)'*y;
    Bv(:,nt,k)=Pv(:,:,k)'*y;
    Be(:,nt,k)=Pomigadot(:,:,k)'*y+Pomiga(:,:,k)'*ydot;
    Ba(:,nt,k)=Pvdot(:,:,k)'*y+Pv(:,:,k)'*ydot;
      
    % ͨ������任���㾲ֹ����ϵ�е����Ĺ��Ծ�
    IT(:,:,k)=SOK(:,:,k)'*IT_Local(:,:,k)*SOK(:,:,k);
 
    % ����ntʱ��,ϵͳ�Ķ���(���ĵ�ƽ��+�������ƽ������ϵ��ת��)
    E_Kin_k(k,nt)=Mass(k)*sum(Bv(:,nt,k).^2)/2+...
        diag(IT_Local(:,:,k)).'*Bw(:,nt,k).^2/2;
     
    % ����ntʱ��,ϵͳ������(����-����)����ע�ⲻͬ
    E_Gra_k(k,nt)=(-Mass(k)*g-FB(3,k))*(pG(3,k)+h_Pre+L_Aux+R_Cha);  
    
end
E_Kin(nt)=sum(E_Kin_k(:,nt));
E_Gra(nt)=sum(E_Gra_k(:,nt));
% ����ֹͣ��־λ
if Bv(1,nt,N)<=1.0e-6
    BvSFlag=1;
end

%% �����ɻ����˶���ѧ����ϵ����a_kl��f_1��
%% �������������,���������,��������ֵ�����֤(����ntʱ�̵�������nt+1ʱ���ü�ydottemp)
% ����������F_Act<6N*1>Ϊ������
F_Act(1:Free,1)=0;  
nt=size(Br,2);
% ��ΪС�Ƕ�ʱsin(theta)=theta???
% =========================================================================  ===================== ���Ľ� ==========================
base_dis=[Bs(:,nt,N);sin(Br(:,nt,N))];
base_acc=[Ba(:,:,N)',Be(:,:,N)']; 
% =========================================================================  ===================== ���Ľ� ==========================
% FCΪˮ��������
vS(1,nt)=Bv(1,nt,N);           vW(1,nt)=1.261;                                                 
vSW(1,nt)=vS(1,nt)-vW(1,nt);    
Re=abs(vSW(1,nt))*LenN(1)/upsilon;
Cf=0.075/(log(Re)-2)^2;
Csf=(105*(ks/LenN(1))^(1/3)-0.64)*1.0e-3; 
FC(1,nt,N)=-sign(vSW(1,nt))*(Cf+Csf)*rou_W*S_Wet*vSW(1,nt).^2/2; 
FC(2,nt,N)=-sign(vSW(2,nt))*rou_W*Cy*Ay*vSW(2,nt).^2/2; 
% FRΪ��ˮ�ָ�������
FR(:,nt)=FR_Factor.*base_dis;
E_FR(nt)=-sum(FR_Factor.*[Bs(:,nt,N).^2;Br(:,nt,N).^2])/2;
F_Act=F_Act+Pv(:,:,N)*FR(1:3,nt)+Pomiga(:,:,N)*FR(4:6,nt)...
    +Pv(:,:,N)*FC(1:3,nt,N);

% =========================================================================  ==================== ���Ľ� ==========================
% ȡ�����Եػ���������ê�����ã���������Լ������
for k=2:N-2
    if  Bv(3,nt,k)<0&&Bs(3,nt,k)<-h_Pre-L_Aux
        if k<=N_Cha+1
            FSC(3,nt,k)=-c_Seabed*(Bs(3,nt,k)+h_Pre+L_Aux)*prod(Len_Cha(1:2));
        else
            FSC(3,nt,k)=-c_Seabed*(Bs(3,nt,k)+h_Pre+L_Aux)*prod(Len_Pre(1:2));
        end
    end
end

% =========================================================================  ==================== ���Ľ� ==========================
FDF(1,1)=-mu*(-Mass(1)*g-FB(3,1));  
if nt>=2
    % F_Hin_KaneΪ��ֵʱ��ʾ����
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
               +Pomiga(:,:,k)*(M(:,k)+MSC(1:3,nt,k)+SOK(:,:,k)*M_Local(:,k)); %�����ֲ�����������
end

% ntʱ��,����������F_Act���������y�ĳ˻�(��ӦԪ����˺����).nt=1:Step_n+1
% ����ĳЩ������ϵͳ����Ҫ�ӳ��������
dK_nt(nt)=Delta_t*F_Act.'*y;
E_Flyl(nt)=sum(dK_nt(1:nt))+E0_Kin;
% a_lpΪ�������������������ʵ���(���ٶ�)�йص�ϵ������,��a_lp
% h_lΪ���������������������йص���,��h_l
a_lp(1:Free,1:Free)=0;  h_l(1:Free,1)=0;     %�Գ�������
for k=1:N
    a_lp=a_lp+Mass(k)*Pv(:,:,k)*Pv(:,:,k)'...
        +Pomiga(:,:,k)*IT(:,:,k)*Pomiga(:,:,k)';
    h_l=h_l+Mass(k)*Pv(:,:,k)*Pvdot(:,:,k)'*y...
        +Pomiga(:,:,k)*IT(:,:,k)*Pomigadot(:,:,k)'*y...
        +Pomiga(:,:,k)*Dual_Mat(Pomiga(:,:,k)'*y)*IT(:,:,k)*(Pomiga(:,:,k)'*y);
end

%% �ﹹ��Լ�����̡�
% ��άê������ϵͳ,��Լ�����̿�дΪB_Sub*y(Var)=g_Ext��Bsub*Ysub=g_Ext;Ysub=y(Var)
% g_Ext_xz��ʾx��z������ⲿԼ������ʾz������ⲿԼ��
g_Ext_xz=zeros(2,1);      g_Ext_z=zeros(1,1); 
g_Extdot_xz=zeros(2,1);   g_Extdot_z=zeros(1,1); 
% ����ⲿԼ������B_Ext���䵼��Bdot
B_Ext1=Pomiga(:,:,1).';  B_Ext2=Pv(:,:,1).';  B_Ext=[B_Ext1; B_Ext2];
% ����������Լ������Bsub[m*r];  r=length(Var)
Bsub_xz=B_Ext([4 6],Var);   Bsub_z=B_Ext(6,Var);
% ���Bsub����������Csub[r*(r-m)];
% Method 1 - Csub is a matrix composed of r-m eigenvectors which ...
% correspond to the zero eigenvalues of Bsub.'*Bsub.
% ����1-Csub����(r-m)��Bsub.'*Bsub������������Ӧ������ʸ��Ϊ�еľ���
[m_xz,r_xz]=size(Bsub_xz);  [C1_xz,D_xz]=eig(Bsub_xz.'*Bsub_xz);      Csub_xz=[];
% ѡ��Bsub.'*Bsub����������Ӧ��������������ɾ���C
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
% ���Bsub�ĵ���Bsubdot
B_Ext1dot=Pomigadot(:,:,1).';     B_Ext2dot=Pvdot(:,:,1).';    
B_Extdot=[B_Ext1dot; B_Ext2dot]; 
Bsubdot_xz=B_Extdot([4 6],Var);
Bsubdot_z=B_Extdot(6,Var);
% ϵͳԼ������,�����󵼺�ΪBsub*ydot(Var)=g_Extdot-Bsubdot*y(Var)

%% �鼯�˶����̣��������ɶ�
% FsumΪ����������йصĹ�����
% AsubΪ�������ɶȺ��ϵͳ�����˶�΢�ַ��̵����ϵ������
% inv(Asub)*FsubΪ������΢�ַ��̵��ұ߲���
% ��������˶���ѧ����ΪAsub*ydot(Var)=Fsub+Bsub.'*lamda...
% ����ͬ����Csub.'�ɵ�Csub.'*Asub*ydot(Var)=Csub.'*Fsub
Fsum=F_Act-h_l;
Asub=a_lp(Var,Var);
Fsub=Fsum(Var);
% ϵͳԼ���������������ϵͳ�˶���ѧ�����������ΪUA*ydot(Var)=UB
UA_xz=[Csub_xz.'*Asub;  Bsub_xz];
UB_xz=[Csub_xz.'*Fsub;  g_Extdot_xz-Bsubdot_xz*y(Var)];
UA_z=[Csub_z.'*Asub;  Bsub_z];
UB_z=[Csub_z.'*Fsub;  g_Extdot_z-Bsubdot_z*y(Var)];
% ��Ϊ��һʱ�̴������ɶȵĹ������ʵĵ���
% F_Hin_KaneΪ��ʱ��ʾ��������ʱ��ʾѹ��
if (nt>=2&&-F_Hin_Kane(1,nt-1,2)>-FDF(1,nt))||Bv(1,nt,1)>1.0e-6
    % ��Ϊ��һʱ�̴������ɶȵĹ������ʵĵ���
    ydottemp=UA_z\UB_z;
    gt_z(:,nt)=Bsub_z*y(Var);
    gtdot_z(:,nt)=Bsub_z*ydottemp+Bsubdot_z*y(Var);
    Error_e1_z(:,nt)=Csub_z.'*Asub*ydottemp-Csub_z.'*Fsub;
    Error_e2_z(:,nt)=Bsub_z*ydottemp-g_Extdot_z+Bsubdot_z*y(Var);
    Error_e_z=[Error_e1_z;  Error_e2_z];
    % ���Լ������Լ�����ص������lamda
    lamda_z(:,nt)=Bsub_xz.'\(Asub*ydottemp-Fsub);
else
    % �÷�����Ҫ��δ���������������������������������������������������������������
    ydottemp=UA_xz\UB_xz;
    % gtΪ��¼��ͬʱ�̵�Լ���� yzΪ����Bsub*ydottemp-g_Extdot+Bsubdot*y(Var)�����
    gt_xz(:,nt)=Bsub_xz*y(Var);
    gtdot_xz(:,nt)=Bsub_xz*ydottemp+Bsubdot_xz*y(Var);
    Error_e1_xz(:,nt)=Csub_xz.'*Asub*ydottemp-Csub_xz.'*Fsub;
    Error_e2_xz(:,nt)=Bsub_xz*ydottemp-g_Extdot_xz+Bsubdot_xz*y(Var);
    Error_e_xz=[Error_e1_xz;  Error_e2_xz];
    % ���Լ������Լ�����ص������lamda
    lamda_xz(:,nt)=Bsub_xz.'\(Asub*ydottemp-Fsub);
end

%% ����Kane�������ӵ�Լ����
% F_conΪ����Լ�������󣬾�����ر任�����ø��ӵ㴦��Լ����
Arem=a_lp(Var_Know,Var);
Frem=Fsum(Var_Know);
% F_ConΪ����Լ����������
F_Con(1:Free,nt)=0; 
F_Con(Var_Know,nt)=Arem*ydottemp-Frem;
% Mcon_3N��Fcon_3N�ֱ�Ϊ����ڸ��������Ĵ���ʩ�ӵ�Լ�����غ�Լ�������ǹ�������
% �����任������������ӵ㴦��Լ�������غ�Լ����
% Pomiga_Reshape*Mcon_3N+Pv_Reshape*Fcon_3N=F_Con(��СΪ6N*1); 
Pomiga_Reshape(:,:)=reshape(Pomiga,6*N,3*N);
Pv_Reshape(:,:)=reshape(Pv,6*N,3*N);  
P_6N=[Pomiga_Reshape Pv_Reshape];
Fcon_6N(:,nt)=P_6N\F_Con(:,nt);
Mcon_3N(:,nt)=Fcon_6N(1:3*N,nt);  Fcon_3N(:,nt)=Fcon_6N(3*N+1:6*N,nt);
% �˶༶��ϵͳ�У�Mcon_k����M_Hin
Mcon_k(:,nt,:)=reshape(Mcon_3N(:,nt),3,N);
Fcon_k(:,nt,:)=reshape(Fcon_3N(:,nt),3,N);
for k=1:N
    F_Hin_Kane(:,nt,k)=sum(Fcon_k(:,nt,k:N),3); 
end

% ��֤����Լ������Լ������
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

