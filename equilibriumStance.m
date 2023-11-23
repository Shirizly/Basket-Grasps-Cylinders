function p3 = equilibriumStance(p1,p2,ph1,ph2,t,hcm)
% this function returns the position of the third finger given two fingers
% position and relative angles on their circular cross-section of the
% object. This doesn't require the grasp to be symmetrical, though the two
% fingers do need to be supporting the base at the same circular
% cross-section (h_12=h_1=h_2).
n1 = [-sin(ph1),-cos(ph1)*cos(t),-cos(ph1)*sin(t)].';
n2 = [-sin(ph2),-cos(ph2)*cos(t),-cos(ph2)*sin(t)].';
n3 = [0,-sin(t),cos(t)].';
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
com = Rt*[0;0;hcm];
N = [n1 n2 n3];
g = [0 0 -1].';
L = -N\g;
A = [n1(1:2) -n2(1:2)];
dp = p2-p1;
sol = A\dp(1:2);
inter1 = p1+sol(1)*n1;
inter2 = p2+sol(2)*n2;
%     h1 = inter1(3)-com(3);
%     h2 = inter2(3)-com(3);
h1 = inter1(3)-com(3);
h2 = inter2(3)-com(3);
lxy = inter1(1:2)-com(1:2);
h3 = (lxy(2)-h1*L(1)*n1(2)-h2*L(2)*n2(2))/(L(3)*n3(2));
inter3 = [lxy;h3]+com;
dist3 = n3.'*inter3;
p3 = double(inter3-n3*dist3);
end