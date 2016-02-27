function y= Dual_Mat(a)
% 该函数为求解(角速度)的对偶矩阵
y=[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];

