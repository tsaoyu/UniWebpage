
% CatenaryShape_HorizontalPCΪ���������������ƽ����״�ĳ���
% ��֪��������L_Pre�������϶������׾���h_Pre������ˮƽͶӰx_Pre������ˮ�е�λ����w_PUW�����������������
% ע������ϵ�����⣺�����������ϵ��������DynRes_AncCha_Anc1�еĹ���ϵ��ͬ
% ������������ϵԭ��λ�����¶����������������������״�����ȣ�
% ����ı任���������н��С���

global  L_Pre  h_Pre  x_Pre  x_PreS  w_PUW  N_Pre  X_HoP  Z_HoP  Tv_Pre  theta_Pre

%% ����fsolve�������������������״�����������ȣ���
% �Ȱ�����ϵԭ��λ�����¶�ê�������״������������ 
% b0Ϊb�ĳ�ֵ
b0=10;       
b=fsolve(@(b)(L_Pre-x_Pre)*(cosh(b)-1)-h_Pre*(sinh(b)-b),b0);
% �������������� 
% a=Th_Pre_Pre/w_PUW, b=x_PreS/a��Ϊ�����������õĲ���; 
Th_Pre=w_PUW*(L_Pre-x_Pre)/(sinh(b)-b); 
Tv_Pre=Th_Pre*sinh(b); 
a=Th_Pre/w_PUW;  
Sc=sinh(b)*a;
x_PreS=a*b;

%% ����ê�����½ӵ��x,z���� 
%  X_HoP, Z_HoP�ֱ�Ϊê�����������x,z�������(������)
X_HoP=zeros(N_Pre+1,1);  Z_HoP=zeros(N_Pre+1,1);
for k=1:N_Pre+1
    % ͨ��s=0��L_Pre�����ê�������½ӵ������
    % ��ע�⣺ê�����εı�Ŵ����¶˿�ʼ1,2,...��
    s=(k-1)*L_Pre/N_Pre;
    if s<=x_Pre-x_PreS
        X_HoP(k)=s;           Z_HoP(k)=0;
    else
        Sc_s=s-(L_Pre-Sc);            xS_s=a*asinh(Sc_s/a);
        X_HoP(k)=x_Pre-(x_PreS-xS_s);        Z_HoP(k)=a*(cosh(xS_s/a)-1);
    end
end
% �����������ʼ�Ƕ�(������)
theta_Pre=zeros(N_Pre,1);
for k=1:N_Pre
    theta_Pre(k)=atan((Z_HoP(k+1)-Z_HoP(k))/(X_HoP(k+1)-X_HoP(k)));
    if abs(theta_Pre(k)-pi/2)<1.0e-6||abs(theta_Pre(k)+pi/2)<1.0e-6
        pause;   disp('��������������ֱ�Σ�')
    end
    %     % ת��Ϊ��x������Ķ۽�(��ֵ)
end
% set(0,'DefaultFigureColor','w');
% set(0,'DefaultAxesFontname','Times New Roman');
% figure;  plot(X_HoP,Z_HoP,'-o','LineWidth',1); 
% title('Exist Horizontal Prevention Cable');