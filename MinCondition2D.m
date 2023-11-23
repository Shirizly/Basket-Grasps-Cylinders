function [flag,val] = MinCondition2D(hcm,r,phi,t1,d,p1,p2,p3)
t = t1;p = phi;
val = (-2*d - cos(t)^2 + r*cos(p)*cos(t)*sin(t))/cos(t);
val2 = - d*cos(t) - (sin(t)*(2*d*sin(t) + r*cos(p)*cos(t)))/cos(t) - (cos(t)*(d + 1))/(d^2*tan(t)^2 + 1)^(1/2);

%%
% com = [-sin(t1);cos(t1)];
% p12 = (p1+p2)/2;
% p12 = p12(2:3);
% p12 = p12-com;
% p3 = p3(2:3);
% p3 = p3-com;
% n3 = [-sin(t1),cos(t1)].';
% n12 = [cos(t1),sin(t1)].';
% N = [n12 n3];
% h = [n12 -n3]\(p3-p12); %distance fingers-to-norm.inters.
% T = p12 + h(1)*n12; %normal intersection
% hcm = -T(2);
% J =  [0 -1;1 0];
% val= hcm - (p12-p3).'*[n3 J*n3]*n12/det(N);
flag = sign(-val)/2+0.5;
% if abs(val-val2)>1E-9
%     disp('error')
% end
end