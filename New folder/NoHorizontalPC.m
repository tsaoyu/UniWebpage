function  F = NoHorizontalPC(x)

global x_Pre  h_Pre  L_Pre
x_PreS=x_Pre;  
%NOHORIZONTALCHAIN Summary of this function goes here
%   Detailed explanation goes here
a=x(1);
theta0=x(2);
F = [a*cosh(x_PreS/a+log(tan(theta0)+sec(theta0)))-a*sec(theta0)-h_Pre;
     a*sinh(x_PreS/a+log(tan(theta0)+sec(theta0)))-a*tan(theta0)-L_Pre];

end

