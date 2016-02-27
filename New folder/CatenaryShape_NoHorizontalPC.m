
% CatenaryShape_NoHorizontalPCΪ���������������ƽ����״�ĳ���
% ��֪������L_Pre�������϶������׾���h_Cha������ˮƽͶӰx_Pre��ê��ˮ�е�λ����w_PUW�����������������
% ע������ϵ�����⣺�����������ϵ��������DynRes_AncCha_Anc1�еĹ���ϵ��ͬ
% ������������ϵԭ��λ�����¶����������������������״�����ȣ�
% ����ı任���������н��С���

global  L_Pre  h_Pre  x_Pre  x_PreS  w_PUW  N_Pre  X_HoP  Z_HoP  Tv_Pre  theta_Pre

%% ����fsolve����������������Է�����(1)(2),solve�����ҵ���ʽ�ⷨ
% (1) a*cosh(x/a+log(tan(theta0)+sec(theta0)))-a*sec(theta0)-z
% (2) a*sinh(x/a+log(tan(theta0)+sec(theta0)))-a*tan(theta0)-L_Pre
% (3) Th*cosh(x/a+log(tan(theta0)+sec(theta0)))-T
% ��xS���¸�ֵ
x_PreS=x_Pre; 
% ��ע������ֵ��ѡ��ʮ����Ҫ�������ʵĳ�ֵ���ܻ�������
sat0=[10, 0];
[sat,fval,exitflag]=fsolve(@ NoHorizontalPC,sat0);
a=sat(1);  theta0=sat(2);   
Th_Pre=a*w_PUW;
T=Th_Pre*cosh(x_PreS/a+log(tan(theta0)+sec(theta0)));
Tv_Pre=sqrt(T^2-Th_Pre^2);

%% �����������½ӵ��x,z���� 
%  X_HoP, Z_HoP�ֱ�Ϊ�������������x,z�������(������)
X_HoP=zeros(N_Pre+1,1);  Z_HoP=zeros(N_Pre+1,1);
for k=1:N_Pre+1
    % ͨ��s=0��L_Pre�����ê�������½ӵ������
    % ��ע�⣺ê�����εı�Ŵ����¶˿�ʼ1,2,...��
    s=(k-1)*L_Pre/N_Pre;
    x_int=fsolve(@(x)a*sinh(x/a+log(tan(theta0)+sec(theta0)))-a*tan(theta0)-s,x_PreS);
    z_int=a*cosh(x_int/a+log(tan(theta0)+sec(theta0)))-a*sec(theta0);
    X_HoP(k)=x_int;
    Z_HoP(k)=z_int;
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
% title('No Horizontal Prevention Cable');  
