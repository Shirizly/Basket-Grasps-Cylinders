syms theta s1 s2 r hcm
T = [cos(theta) -sin(theta);sin(theta) cos(theta)];
% v1 = sym('v1',[2,1]);
% v2 = sym('v2',[2,1]);
% t1 = sym('t1',[2,1]);
% t2 = sym('t2',[2,1]);
v1 = [r;-hcm];
v2 = [-r;-hcm];
t1 = [-1;0];
t2 = [0;1];
syms p_c11 p_c12 p_c21 p_c22

p_c1 = [p_c11;p_c12];
p_c2 = [p_c21;p_c22];

exp1 = p_c2 - T*(v1+s1*t1);
exp2 = p_c1 - T*(v2+s2*t2);

sol = solve(exp1==exp2,[s1,s2],'ReturnConditions',true);

s1f = simplify(sol.s1)
s2f = simplify(sol.s2)
%%
p_rel = simplify(p_c2 - T*(v1+s1f*t1))
% dpdtheta = simplify(diff(p_rel(2),theta));
% dp2dtheta2 = simplify(diff(dpdtheta,theta));
dpdtheta = p_c11*cos(2*theta) - p_c21*cos(2*theta) + p_c12*sin(2*theta) - p_c22*sin(2*theta) + r*cos(theta) - hcm*sin(theta);
d2pdtheta2 = 2*p_c12*cos(2*theta) - 2*p_c22*cos(2*theta) - 2*p_c11*sin(2*theta) + 2*p_c21*sin(2*theta) - hcm*cos(theta) - r*sin(theta);

sol = solve(dpdtheta,theta,'ReturnConditions',true)%,'MaxDegree',4)
%p_c21*z^4*1i - p_c11*z^4*1i - p_c12*z^4 + p_c22*z^4 - r*z^3*1i + hcm*z^3 - r*z*1i - hcm*z + p_c21*1i - p_c22 - p_c11*1i + p_c12
%p_c21*z^4*1i - p_c11*z^4*1i - p_c12*z^4 + p_c22*z^4 - r*z^3*1i + hcm*z^3 - r*z*1i - hcm*z + p_c21*1i - p_c22 - p_c11*1i + p_c12
polyRep = [p_c21*1i - p_c11*1i - p_c12 + p_c22 , - r*1i + hcm , 0 , - r*1i - hcm , p_c21*1i - p_c22 - p_c11*1i + p_c12];

%%
matlabFunction(s1f,s2f,p_rel,dpdtheta,dp2dtheta2,'File','DSStateFromTheta')

%%
matlabFunction(sol.theta,'File','DSRecEscape')