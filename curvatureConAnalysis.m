syms t p d r hcm
% hcm = 1;

a = [0,-sin(t),cos(t)].';
% contact normals
n1 = -[sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n2 = -[-sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n3 = a;

%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = simplify(-N\g);

%% points of contact, frame origin is at intersection of l1,l2
ro1 = -n1*r;
ro2 = -n2*r;
rocm = -n3*d;
n12 = [0;cos(t);sin(t)]; %projection of the normals 1,2 to the cross section that includes p3, n3, and COM
rob = n12*d*sin(t)/cos(t); %intersection of n3 with the plane that includes n1 and n2
robn = sign(d)*[0;cos(t);sin(t)];
normrob = d*tan(t);
ro3 = rob-n3*(hcm+d);

%% testing point for expressions
% cross1 = sTra(-CrMat(ro1)*CrMat(n1))*l3,20)
% cross2 = sTra(-CrMat(ro2)*CrMat(n2))*l3,20)
% cross1+cross2)

%% curvature form
sandwich_mat_1 = [eye(3) -CrMat(ro1);zeros(3) CrMat(n1)];
sandwich_mat_2 = [eye(3) -CrMat(ro2);zeros(3) CrMat(n2)];
sandwich_mat_3 = [eye(3) -CrMat(ro3);zeros(3) CrMat(n3)];

U1 = [zeros(3,1);a];
l2 = robn;%[0;cos(t);sin(t)];
U2 = [zeros(3,1);l2];
x = [1;0;0];
l3 = x;
U3 = [cross(rob,l3);l3];
e = [0;0;1];
U = [U2,U3];
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
pb3 = CrMat(ro3).'*CrMat(n3);
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
% kright13 = [n3*(-d*tan(t)+cos(p)*r);n3*cos(p)];
% 
% kright2 = sandwich_mat_2*U;
% kright21 = [cross(n2,a)*r;cross(n2,a)];
% kright22 = [cross(n2,l2)*r;cross(n2,l2)];
% kright23 = [n3*(-d*tan(t)+cos(p)*r);a*cos(p)];
% % kright2-[kright21 kright22 kright23];
% 
% kright3 = sandwich_mat_3*U;
% kright31 = [-cross(ro3,a);cross(a,a)];
% kright32 = [-cross(ro3,robn);cross(a,l2)];
% kright33 = [(d+hcm)*robn;robn];
% % kright3-[kright31 kright32 kright33];



%%
t1 = cross(n1,a);t2 = cross(n2,a);
% curv1sandU = [zeros(3,1) sin(p)*a -cos(p)*a; -r*t1 -r*sin(p)*a (normrob-r*cos(p))*a];
% curv2sandU = [zeros(3,1) -sin(p)*a -cos(p)*a; -r*t2 r*sin(p)*a (normrob-r*cos(p))*a];
% curv3sandU = [zeros(3,1) l3 -robn; normrob*l3 (d+hcm)*l3 (d+hcm)*robn];


curv1sandU = [sin(p)*a -cos(p)*a;  -r*sin(p)*a (normrob+r*cos(p))*a];
curv2sandU = [-sin(p)*a -cos(p)*a;  r*sin(p)*a (normrob+r*cos(p))*a];
curv3sandU = [l3 -robn;  (d+hcm)*l3 (d+hcm)*robn];

% curv1sandU = [-cos(p)*a;(normrob+r*cos(p))*a];
% curv2sandU = [ -cos(p)*a; (normrob+r*cos(p))*a];
% curv3sandU = [-robn;(d+hcm)*robn];

sum1 = sandwich_mat_1.'*curv1sandU;
sum2 = sandwich_mat_2.'*curv2sandU;
sum3 = sandwich_mat_3.'*curv3sandU;

c1 = 1;
c2 = 1;
c3 = sqrt(1+((d*sin(t))/cos(t))^2);
G =  L(1)/c1*(sum1-P1U)+L(2)/c2*(sum2-P2U)+L(3)/c3*(sum3-P3U)-D2U*U;
G = U.'*G;
pretty(simplify(G))
%% according to Elon's development
u3 = [1;0;0];
k = [t;p;d];
k(1) = u3.'*CrMat(2*rob-r*n1)*CrMat(n1)*u3;
k(2) = u3.'*CrMat(2*rob-r*n2)*CrMat(n2)*u3;
d3 = d+hcm;
crau3 = cross(a,u3)
c3 = sqrt(1+d^2*tan(t)^2);
k(3) = -1/c3*d3*crau3.'*crau3;
k = simplify(k)
Ku = L(1)*k(1) + L(2)*k(2)+L(3)*k(3)-U3.'*D2U*U3;
Ku = simplify(Ku)
%% according to Elon, new
k1 = [t,t;t,t];
w = [l2,l3];
U = [U2 U3];
for i=1:2
    k1(i,i) = w(:,i).'*CrMat(-2*rob+ro1).'*CrMat(n1)*w(:,i);
end
% for i=1:2
%     k1(i,3-i) = l2.'*CrMat(-2*rob+ro1).'*CrMat(n1)*l3;
% end
k1(1,2) = l2.'*CrMat(-2*rob+ro1).'*CrMat(n1)*l3;
k1(2,1) = l3.'*CrMat(-2*rob+ro1).'*CrMat(n1)*l2;

k2 = [t,t;t,t];
for i=1:2
    k2(i,i) = w(:,i).'*CrMat(-2*rob+ro2).'*CrMat(n2)*w(:,i);
end
% for i=1:2
%     k2(i,3-i) = l2.'*CrMat(-2*rob+ro2).'*CrMat(n2)*l3;
% end
k2(1,2) = l2.'*CrMat(-2*rob+ro2).'*CrMat(n2)*l3;
k2(2,1) = l3.'*CrMat(-2*rob+ro2).'*CrMat(n2)*l2;

k3 = [t,t;t,t];
for i=1:2
    k3(i,i) = -(d+hcm);
end
for i=1:2
    k3(i,3-i) = 0;
end
k3 = k3/sqrt(1+normrob^2);
k12 =  simplify(k1 + k2);
Ku12 =  simplify(L(1)*k1 + L(2)*k2);
Ku3 = simplify(L(3)*k3);
Kuu = simplify(U.'*D2U*U);
Ku = Ku12 + Ku3 - Kuu;
Ku = simplify(Ku)

%% redevelop
% P1 = CrMat(ro1).'*CrMat(n1);
% P1 = sTra(CrMat(ro1).'*CrMat(n1));
% WP1W1 = -w.'*P1*w;WP1W1 = simplify(WP1W1)
WP1W = [r*cross(n1,robn).'*cross(n1,robn) r*cross(n1,robn).'*cross(n1,x);...
    r*cross(n1,robn).'*cross(n1,x) r*cross(n1,x).'*cross(n1,x)];
WP1W = simplify(WP1W);
WS1W = -[2*r*cross(n1,robn).'*cross(n1,robn) cross(n1,robn).'*cross(rob-2*ro1,x);...
    cross(n1,robn).'*cross(rob-2*ro1,x) 2*cross(rob-ro1,x).'*cross(n1,x)];
k1 = simplify(WS1W+WP1W);

k1t = [-r*cross(n1,robn).'*cross(n1,robn) cross(n1,robn).'*cross(ro1-rob,x);...
    cross(n1,robn).'*cross(ro1-rob,x) cross(ro1-2*rob,x).'*cross(n1,x)];
simplify(k1t-k1)


P2 = sTra(CrMat(ro2).'*CrMat(n2));
% WP2Wt = -w.'*P2*w;WP2Wt = simplify(WP2Wt)
WP2W = [r*cross(n2,robn).'*cross(n2,robn) r*cross(n2,robn).'*cross(n2,x);...
    r*cross(n2,robn).'*cross(n2,x) r*cross(n2,x).'*cross(n2,x)];
WP2W = simplify(WP2W);
WS2W = -[2*r*cross(n2,robn).'*cross(n2,robn) cross(n2,robn).'*cross(rob-2*ro2,x);...
    cross(n2,robn).'*cross(rob-2*ro2,x) 2*cross(rob-ro2,x).'*cross(n2,x)];
k2 = simplify(WS2W+WP2W);
k12 = simplify(k1+k2);
k12t = simplify([-2*r*cross(n1,robn).'*cross(n1,robn) 0;0 cross(ro1+ro2-4*rob,x).'*cross(n1,x)]);
simplify(k12t-k12)

sandwich_mat_3 = [eye(3) -CrMat(ro3);zeros(3) CrMat(n3)];
WS3W = U.'*sandwich_mat_3.'*[zeros(3) -eye(3);-eye(3) zeros(3)]*sandwich_mat_3*U;
WS3W = simplify(WS3W);
P3 = sTra(CrMat(ro3).'*CrMat(n3));
WP3W = -w.'*P3*w;WP3W = simplify(WP3W);
k3 = WS3W + WP3W;
k3t = -(d+hcm)*eye(2);
k3 = k3/sqrt(1+normrob^2);

G = Kuu-L(1)*k12t - L(3)*k3;
G = subs(G,sign(d)^2,1);
G = simplify(G,1000)
 
%%

k12 = [-2*r*cross(n1,robn).'*cross(n1,robn) 0;0 cross(ro1+ro2-4*rob,x).'*cross(n1,x)];

k3 = -(d+hcm)*eye(2);


k3 = k3/sqrt(1+normrob^2);

Ku12 =  simplify(k12);
Ku3 = simplify(k3);
% Kuu = -U.'*D2U*U;
Kuu = d*cos(t)*eye(2);
Ku = simplify(Kuu- L(1)*Ku12 - L(3)*Ku3 );

%%
val = eig(Ku);



function Ax = CrMat(a)
Ax = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end
function As = sTra(A)
As = 0.5*(A+A.');
end