syms t p ds r %hcm
hcm = 1;
a = [0,-sin(t),cos(t)].';
% contact normals
n1 = -[sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n2 = -[-sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n3 = a;

%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = simplify(-N\g)

%% points of contact, frame origin is at intersection of l1,l2
ro1 = -n1*r;
ro2 = -n2*r;
rocm = -n3*ds;
n12 = [0,cos(t),sin(t)].'; %projection of the normals 1,2 to the cross section that includes p3, n3, and COM
rob = n12*ds*sin(t)/cos(t); %intersection of n3 with the plane that includes n1 and n2
robn = [0;cos(t);sin(t)];
normrob = ds*tan(t);
ro3 = rob-n3*(hcm+ds);

%% testing point for expressions
% cross1 = simplify(sTra(-CrMat(ro1)*CrMat(n1))*l3,20)
% cross2 = simplify(sTra(-CrMat(ro2)*CrMat(n2))*l3,20)
% simplify(cross1+cross2)

%% curvature form
sandwich_mat_1 = [eye(3) -CrMat(ro1);zeros(3) CrMat(n1)];
sandwich_mat_2 = [eye(3) -CrMat(ro2);zeros(3) CrMat(n2)];
sandwich_mat_3 = [eye(3) -CrMat(ro3);zeros(3) CrMat(n3)];

U1 = [zeros(3,1);a];
l2 = [0;cos(t);sin(t)];
U2 = [zeros(3,1);l2];
l3 = [1;0;0];
U3 = [cross(rob,l3);l3];
e = [0;0;1];
U = [U1 U2 U3];
%% star matrices
% star = simplify(sTra(CrMat(rocm)*CrMat(e)))
% D2U = [zeros(3) zeros(3);zeros(3) star]
% star1 = simplify([zeros(3) zeros(3);zeros(3) simplify(sTra(CrMat(ro1).'*CrMat(n1)))])
% star2 = simplify([zeros(3) zeros(3);zeros(3) simplify(sTra(CrMat(ro2).'*CrMat(n2)))])
% star3 = simplify([zeros(3) zeros(3);zeros(3) simplify(sTra(CrMat(ro3).'*CrMat(n3)))])
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
% D2U_q = simplify(hessian(rho_cm.'*e,q))
% D2U_0 = subs(D2U_q,th,[0;0;0])
% pretty(D2U_q)
% pretty(D2U_0)
starcm = diag([-rocm(3),-rocm(3),0]);
D2U = [zeros(3,6); zeros(3),starcm]
%% contact matrices
conblock1 = sTra(CrMat(ro1).'*CrMat(n1));
con1 = [zeros(3,6); zeros(3), conblock1]
conblock2 = sTra(CrMat(ro2).'*CrMat(n2));
con2 = [zeros(3,6); zeros(3), conblock2]
conblock3 = sTra(CrMat(ro3).'*CrMat(n3));
con3 = [zeros(3,6); zeros(3), conblock3]
%% version with the symmetrization effect
% conblock1 = CrMat(ro1).'*CrMat(n1);
% con1 = [zeros(3,6); zeros(3), conblock1]
% conblock2 = CrMat(ro2).'*CrMat(n2);
% con2 = [zeros(3,6); zeros(3), conblock2]
% conblock3 = CrMat(ro3).'*CrMat(n3);
% con3 = [zeros(3,6); zeros(3), conblock3]
%%
conU1 = simplify(con1*U)
conU2 = simplify(con2*U)
conU3 = simplify(con3*U)
%% middle steps in computation
kright1 = simplify(sandwich_mat_1*U)

kright11 = simplify([cross(n1,n3)*r;cross(n1,n3)])
kright12 = simplify([cross(n1,l2)*r;cross(n1,l2)])
kright13 = simplify([n3*(-ds*tan(t)+cos(p)*r);n3*cos(p)])

kright2 = simplify(sandwich_mat_2*U)
kright21 = simplify([cross(n2,a)*r;cross(n2,a)])
kright22 = simplify([cross(n2,l2)*r;cross(n2,l2)])
kright23 = simplify([n3*(-ds*tan(t)+cos(p)*r);a*cos(p)])
simplify(kright2-[kright21 kright22 kright23])

kright3 = simplify(sandwich_mat_3*U)
kright31 = simplify([-cross(ro3,a);cross(a,a)])
kright32 = simplify([-cross(ro3,robn);cross(a,l2)])
kright33 = [(ds+hcm)*robn;robn]
simplify(kright3-[kright31 kright32 kright33])



%%
t1 = cross(n1,a);t2 = cross(n2,a);
curve1 = [zeros(3,1) sin(p)*a -cos(p)*a; -r*t1 -r*sin(p)*a (normrob+r*cos(p))*a];
curve2 = [zeros(3,1) -sin(p)*a -cos(p)*a; -r*t2 r*sin(p)*a (normrob+r*cos(p))*a];
curve3 = [zeros(3,1) l3 -robn; normrob*l3 (ds+hcm)*l3 (ds+hcm)*robn];

sum1 = simplify(sandwich_mat_1.'*curve1)
sum2 = simplify(sandwich_mat_2.'*curve2)
sum3 = simplify(sandwich_mat_3.'*curve3)

c1 = 1;
c2 = 1;
c3 = sqrt(1+((ds*sin(t))/cos(t))^2);
G =  simplify(D2U*U+L(1)/c1*(sum1-conU1)+L(2)/c2*(sum2-conU2)+L(3)/c3*(sum3-conU3));
G = simplify(U.'*G)


function Ax = CrMat(a)
Ax = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end
function As = sTra(A)
As = 0.5*(A+A.');
end