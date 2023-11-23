syms t1 hcm ds r phi
hcm = 1
com = hcm*[-sin(t1);cos(t1)];
p12 = com*(hcm+ds)/hcm+r*(-cos(phi))*[-cos(t1);-sin(t1)];
n12 = [cos(t1);sin(t1)];
n3 = [-sin(t1);cos(t1)];
I2 = p12(2)+(p12(1)-com(1))/(-n12(1))*n12(2);
I = [com(1);I2];
hp = I2-com(2);
p3 = simplify([com(1);I2]-n3*(hcm+ds))
gravity = [0;-1];
J = [0 1;-1 0];
grmat = [n12 n3 gravity;
         -(com-p12).'*J*n12 -(com-p3).'*J*n3 0];
condEq = simplify(det(grmat)); % equilibrium condition, fulfilled automatically by construction
N = [n12 n3];
condMin = simplify(hp - (p12-p3).'*[n3 J*n3]*n12)
pretty(condMin)


%%
p1  = [-2;1];
p2 = [2;1];
syms x y

%%
t1 = pi()/4;
hcm = 1;
h = 2*hcm;
r = 0.5;
ds = 0.3;
phi = pi()*3/4;
p12n = subs(p12);
p3n = subs(p3);
comn = subs(com);
In = subs(I);
n12n = subs(n12);
n3n = subs(n3);
figure
hold on

%plot rectangle view of the cylinder
t = t1;
v1 = -r*[cos(t);sin(t)];
v2 = r*[cos(t);sin(t)];
v3 = r*[cos(t);sin(t)]+h*[-sin(t);cos(t)];
v4 = -r*[cos(t);sin(t)]+h*[-sin(t);cos(t)];
vertex = [v1 v2 v3 v4 v1];
plot(vertex(1,:),vertex(2,:),'b');
% plot axis of cylinder
axisends = [v1+v2 (v3+v4)/2];
plot(axisends(1,:),axisends(2,:),'b--');
% mark com
plot(comn(1),comn(2),'*r');
plot([comn(1),comn(1)],[0,comn(2)*2],'r--');
%plot contact points
plot(p12n(1),p12n(2),'ok','markerSize',4,'lineWidth',4);
plot(p3n(1),p3n(2),'ok','markerSize',4,'lineWidth',4);
% plot force lines
f12 = [p12n p12n+n12n];
plot(f12(1,:),f12(2,:),'k--');

f3 = [p3n p3n+2*n3n];
plot(f3(1,:),f3(2,:),'k--');
xlabel('Y')
ylabel('Z')

plot(In(1),In(2),'*g','markerSize',4,'lineWidth',4);
axis equal
double(subs(condMin))
double(subs(hp))