syms t p ds r 
hcm = 1;
a = [0,-sin(t),cos(t)].';
% contact normals
n1 = -[sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n2 = -[-sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n3 = a;
n12 = [0,cos(t),sin(t)].';
rob = n12*ds*sin(t)/cos(t);
ro3 = rob-n3*(hcm+ds);
ro1 = -n1*r;
ro2 = -n2*r;
com = -n3*ds;
p12 = (ro1+ro2)/2;
p12 = p12-com;
p3 = ro3-com;
p12 = p12(2:3);
p3 = p3(2:3);
n3 = [-sin(t),cos(t)].';
n12 = [cos(t),sin(t)].';
N = [n12 n3];
h = [n12 -n3]\(p3-p12); %distance fingers-to-norm.inters.
T = p12 + h(1)*n12; %normal intersection
hcm = -T(2);
J =  [0 -1;1 0];
val= hcm - (p12-p3).'*[n3 J*n3]*n12/det(N);