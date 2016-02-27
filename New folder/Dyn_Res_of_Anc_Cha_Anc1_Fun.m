function [out1,out2,out3] =Dyn_Res_of_Anc_Cha_Anc1_Fun(t,col,flag)
% ������Ϊê������ϵͳ�����񡪿������ֺ���(�����Ļ����˶�ѧ����+ŷ����������+�ƶ�λ�Ʒ���)

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

if nargin < 3 || isempty(flag)
    
    %% �ɻ��ֱ�������ֵ���¸����󣨹������ʡ�ŷ���������ƶ�λ�ƣ�
    ntemp=size(Var,2);  
    % ����δ֪������������
    y(Var)=col(1:ntemp);  
    % ����ŷ����������
    epsi=reshape(col(ntemp+1:ntemp+N*4),4,N);
    % �������λ��ʸ������
    for i=1:length(Varmove_row)
        s(Varmove_row(i),Varmove_col(i))=col(ntemp+4*N+i);
    end

    %% �������Ա任���󼰾��Ա任�����(��Ա任������ŷ�������й�)
    % ע���ʱ�ļ��㷽ʽ���������м���SJK0��ͬ
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
    
    %% �����ƫ���ٶȡ�ƫ�ٶȼ����ǵĵ������������ʡ����Ա任���󡢾��Ա任����ĵ��������ӵ�ʸ��������ʸ�����ƶ��������ƶ������йأ�
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
        Omiga(:,k)=Pomiga(1:Free,:,k)'*y(1:Free);
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
    
    %% �����ɻ����˶���ѧ����ϵ����a_kl��f_1��
    %% ��������������͹����������(��Dyn_Res_of_Anc_Cha_Anc1_Fun_Bk�м����ظ�)
    % ����������F_Act<6N*1>Ϊ������
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
            +Pomiga(:,:,k)*(M(:,k)+MSC(1:3,nt,k)+SOK(:,:,k)*M_Local(:,k)); %�����ֲ�����������
    end

    % ͨ������任���㾲ֹ����ϵ�е����Ĺ��Ծ�
    for k=1:N
        % Ϊʲô�������ı任��ʽ???����-��p26
        IT(:,:,k)=SOK(:,:,k)'*IT_Local(:,:,k)*SOK(:,:,k);
    end
    % ������������
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
        dy=UA_z\UB_z;
    else
        % �÷�����Ҫ��δ���������������������������������������������������������������
        dy=UA_xz\UB_xz;  
    end
    
    % depsiΪŷ�������任����ĵ���(��ŷ�������͹��������й�)
    for k=1:N
        depsi(4,k)=-epsi(1:3,k)'*y(3*k-2:3*k)/2;
        depsi(1:3,k)=(-Dual_Mat(epsi(1:3,k))+epsi(4,k)*E)*y(3*k-2:3*k)/2;
    end
    ds=y(Varmove);
    % ��������Ĳ���,������������������ĸ���Ϊ(6N+4N+3N)
    out1=[dy; reshape(depsi,4*N,1); ds];
    
else
    
    switch(flag)
        case 'init'
            % ��ʼ������,��ʱ����[tspan,y0,options]=odefun([],[],'init')
            out1=0:Delta_t:Step_n*Delta_t;
            out2=col0;
            out3=odeset('reltol',0.001,'maxstep',Delta_t,'InitialStep',Delta_t/2);
        otherwise
            error(['Unknown flag ''' flag '''.']);
    end
end