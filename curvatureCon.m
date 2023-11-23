function [flag,val] = curvatureCon(t,p,ds,r,hcm)
debug = 0;
if 1
a = [0,-sin(t),cos(t)].';
% contact normals
n1 = -[sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n2 = -[-sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n3 = a;
x = [1;0;0];
%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = -N\g;

%% points of contact, frame origin is at intersection of l1,l2
ro1 = -n1*r;
ro2 = -n2*r;
rocm = -n3*ds;
n12 = [0;cos(t);sin(t)]; %projection of the normals 1,2 to the cross section that includes p3, n3, and COM
rob = n12*ds*sin(t)/cos(t); %intersection of n3 with the plane that includes n1 and n2
robn = [0;cos(t);sin(t)];
normrob = ds*tan(t);
ro3 = rob-n3*(hcm+ds);

%%
if debug
figure
hold on
plot3(ro1(1),ro1(2),ro1(3),'ok')
plot3(ro2(1),ro2(2),ro2(3),'ok')
plot3(ro3(1),ro3(2),ro3(3),'ok')
plot3(rocm(1),rocm(2),rocm(3),'+r')
axisends = [ro1 ro1+r*n1];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k-','lineWidth',2);
axisends = [ro2 ro2+r*n2];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k-','lineWidth',2);
axisends = [ro3 ro3+(hcm+ds)*n3];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k-','lineWidth',2);
axisends = [zeros(3,1) rob];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'g-','lineWidth',2);
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
end
%% testing point for expressions
% cross1 = sTra(-CrMat(ro1)*CrMat(n1))*l3,20)
% cross2 = sTra(-CrMat(ro2)*CrMat(n2))*l3,20)
% cross1+cross2)

%% curvature form
sandwich_mat_1 = [eye(3) -CrMat(ro1);zeros(3) CrMat(n1)];
sandwich_mat_2 = [eye(3) -CrMat(ro2);zeros(3) CrMat(n2)];
sandwich_mat_3 = [eye(3) -CrMat(ro3);zeros(3) CrMat(n3)];

U1 = [zeros(3,1);a];
l2 = [0;cos(t);sin(t)];
U2 = [zeros(3,1);l2];
l3 = x;
U3 = [cross(rob,l3);l3];
e = [0;0;1];
U = [U3];
% sandwich_mat_1 = [eye(3) -CrMat(ro1);zeros(3) CrMat(a)];
% sandwich_mat_2 = [eye(3) -CrMat(ro2);zeros(3) CrMat(l2)];
% sandwich_mat_3 = [eye(3) -CrMat(ro3);zeros(3) CrMat(l3)];
%% star matrix and 2nd derivative
star = sTra(CrMat(rocm)*CrMat(e));
D2U = [zeros(3) zeros(3);zeros(3) star];


%% D2U
% syms a b c x y z xcm ycm zcm
% d = [x;y;z]
% e = [0;0;1];
% rho_cm_0 = [xcm;ycm;zcm]
% pretty(rho_cm_0)
% Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
% pretty(Rx)
% Ry = [cos(b) 0 sin(b);0 1 0; -sin(b) 0 cos(b)]
% pretty(Ry)
% Rz = [cos(c) -sin(c) 0;sin(c) cos(c) 0; 0 0 1]
% pretty(Rz)
% R3a = Rz*Ry*Rx;
% rho_cm = R3a*rho_cm_0+d;
% pretty(rho_cm)
% 
% q = [x;y;z;a;b;c];
% th = [a;b;c];
% U = rho_cm.'*e
% grad = gradient(rho_cm.'*e,q)
% pretty(grad)
% grad_0 = subs(grad,th,[0;0;0])
% pretty(grad_0)
% D2U_q = hessian(rho_cm.'*e,q))
% D2U_0 = subs(D2U_q,th,[0;0;0])
% pretty(D2U_q)
% pretty(D2U_0)
% starcm = diag([-rocm(3),-rocm(3),0]);
% D2U = [zeros(3,6); zeros(3),starcm];
%% contact matrices
Pblock1 = sTra(CrMat(ro1).'*CrMat(n1));
P1 = [zeros(3,6); zeros(3), Pblock1];
Pblock2 = sTra(CrMat(ro2).'*CrMat(n2));
P2 = [zeros(3,6); zeros(3), Pblock2];
Pblock3 = sTra(CrMat(ro3).'*CrMat(n3));
P3 = [zeros(3,6); zeros(3), Pblock3];

P1U = P1*U;
P2U = P2*U;
P3U = P3*U;
%% middle steps in computation
% kright1 = sandwich_mat_1*U;
% 
% kright11 = [cross(n1,n3)*r;cross(n1,n3)];
% kright12 = [cross(n1,l2)*r;cross(n1,l2)];
% kright13 = [n3*(-ds*tan(t)+cos(p)*r);n3*cos(p)];
% 
% kright2 = sandwich_mat_2*U;
% kright21 = [cross(n2,a)*r;cross(n2,a)];
% kright22 = [cross(n2,l2)*r;cross(n2,l2)];
% kright23 = [n3*(-ds*tan(t)+cos(p)*r);a*cos(p)];
% % kright2-[kright21 kright22 kright23];
% 
% kright3 = sandwich_mat_3*U;
% kright31 = [-cross(ro3,a);cross(a,a)];
% kright32 = [-cross(ro3,robn);cross(a,l2)];
% kright33 = [(ds+hcm)*robn;robn];
% % kright3-[kright31 kright32 kright33];



%%
% t1 = cross(n1,a);t2 = cross(n2,a);
% % curv1sandU = [zeros(3,1) sin(p)*a -cos(p)*a; -r*t1 -r*sin(p)*a (normrob+r*cos(p))*a];
% % curv2sandU = [zeros(3,1) -sin(p)*a -cos(p)*a; -r*t2 r*sin(p)*a (normrob+r*cos(p))*a];
% % curv3sandU = [zeros(3,1) l3 -robn; normrob*l3 (ds+hcm)*l3 (ds+hcm)*robn];
% 
% 
% % curv1sandU = [sin(p)*a -cos(p)*a;  -r*sin(p)*a (normrob+r*cos(p))*a];
% % curv2sandU = [-sin(p)*a -cos(p)*a;  r*sin(p)*a (normrob+r*cos(p))*a];
% % curv3sandU = [l3 -robn;  (ds+hcm)*l3 (ds+hcm)*robn];
% 
% curv1sandU = [-cos(p)*a;-(normrob-r*cos(p))*a];
% curv2sandU = [ -cos(p)*a; -(normrob-r*cos(p))*a];
% curv3sandU = [-robn;(ds+hcm)*robn];
% 
% 
% sum1 = sandwich_mat_1.'*curv1sandU;
% sum2 = sandwich_mat_2.'*curv2sandU;
% sum3 = sandwich_mat_3.'*curv3sandU;
% 
% c1 = 1;
% c2 = 1;
% c3 = sqrt(1+((ds*sin(t))/cos(t))^2);
% G =  L(1)/c1*(sum1-P1U)+L(2)/c2*(sum2-P2U)+L(3)/c3*(sum3-P3U)-D2U*U;
% G = U.'*G;

%% Elon's redevelopment

w = [l2,l3];
U = [U2 U3];
% k1 = zeros(2,2);
% for i=1:2
%     k1(i,i) = w(:,i).'*CrMat(-2*rob+ro1).'*CrMat(n1)*w(:,i);
% end
% for i=1:2
%     k1(i,3-i) = l2.'*CrMat(-2*rob+ro1).'*CrMat(n1)*l3;
% end
% k2 = zeros(2,2);
% for i=1:2
%     k2(i,i) = w(:,i).'*CrMat(-2*rob+ro2).'*CrMat(n2)*w(:,i);
% end
% for i=1:2
%     k2(i,3-i) = l2.'*CrMat(-2*rob+ro2).'*CrMat(n2)*l3;
% end
k12 = [-2*r*cross(n1,robn).'*cross(n1,robn) 0;0 cross(ro1+ro2-4*rob,x).'*cross(n1,x)];

k3 = -(ds+hcm)*eye(2);


k3 = k3/sqrt(1+normrob^2);

Ku12 =  L(1)*k12;
Ku3 = L(3)*k3;
% Kuu = -U.'*D2U*U;
Kuu = ds*cos(t)*eye(2);
Ku = Kuu - Ku12 - Ku3;
end


val = eig(Ku);

if all(val>=0)
    flag = 1;
else
    flag = 0;
end
end

function Ax = CrMat(a)
Ax = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end
function As = sTra(A)
As = 0.5*(A+A.');
end