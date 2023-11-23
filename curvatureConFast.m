function [flag,Eval] = curvatureConFast(t,p,d,r,hcm)
%check the stability using the relative curvature form of U.
if d>0 % known to be sufficient condition for stability.
    flag = 1;
    Eval = [];
    return
end
Eval(1) = d*cos(t) + (cos(t)*(d + hcm))/(d^2*tan(t)^2 + 1)^(1/2) - r*sin(p)*tan(p)*sin(t); 
Eval(2) = d*cos(t) + tan(t)*(2*d*sin(t) - r*cos(p)*cos(t)) + (cos(t)*(d + hcm))/(d^2*tan(t)^2 + 1)^(1/2);
flag = 1;
if Eval(2)<0||Eval(1)<0
    flag = 0;
end
end