syms t p d r 
hcm = 1;

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
ro3 = -n3*(hcm+d);
rob = -n12*d*sin(t)/cos(t); %intersection of n1 and n2, ds above the com along a
ro1 = rob-n1*r;
ro2 = rob-n2*r;
rocm = -ds/cos(t)*[0;0;1];
n12 = [0;cos(t);sin(t)]; %projection of the normals 1,2 to the cross section that includes p3, n3, and COM

robn = -[0;cos(t);sin(t)];
normrob = d*tan(t);


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
l3 = [1;0;0];
U3 = [cross(rob,l3);l3];
e = [0;0;1];
U = [U3];

%% star matrix and 2nd derivative
star = sTra(CrMat(rocm)*CrMat(e));
D2U = [zeros(3) zeros(3);zeros(3) star];

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
% curv1sandU = [zeros(3,1) sin(p)*a -cos(p)*a; -r*t1 -r*sin(p)*a (normrob+r*cos(p))*a];
% curv2sandU = [zeros(3,1) -sin(p)*a -cos(p)*a; -r*t2 r*sin(p)*a (normrob+r*cos(p))*a];
% curv3sandU = [zeros(3,1) l3 -robn; normrob*l3 (d+hcm)*l3 (d+hcm)*robn];


% curv1sandU = [sin(p)*a -cos(p)*a;  -r*sin(p)*a (normrob+r*cos(p))*a];
% curv2sandU = [-sin(p)*a -cos(p)*a;  r*sin(p)*a (normrob+r*cos(p))*a];
% curv3sandU = [l3 -robn;  (d+hcm)*l3 (d+hcm)*robn];

curv1sandU = [-cos(p)*a;(normrob+r*cos(p))*a];
curv2sandU = [ -cos(p)*a; (normrob+r*cos(p))*a];
curv3sandU = [-robn;(d+hcm)*robn];

sum1 = sandwich_mat_1.'*curv1sandU;
sum2 = sandwich_mat_2.'*curv2sandU;
sum3 = sandwich_mat_3.'*curv3sandU;

c1 = 1;
c2 = 1;
c3 = sqrt(1+((d*sin(t))/cos(t))^2);
G =  L(1)/c1*(sum1-P1U)+L(2)/c2*(sum2-P2U)+L(3)/c3*(sum3-P3U)-D2U*U;
G = U.'*G;
disp(G)
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
Ku = L(1)*k(1) + L(2)*k(2)+L(3)*k(3)+U3.'*D2U*U3;
%% Redevelop myself
u3 = [1;0;0];
kA = [t;p;d];
% contact matrices results:
U3p1U3 = -r*cos(p)^2;
U3p2U3 = -r*cos(p)^2;
U3p3U3 = -(d+1);

etta1 = [cross(ro1,n1);n1]
etta2 = [cross(ro2,n2);n2]
etta3 = simplify([cross(ro3,n3);n3])

%Curvature forms:
kA(1) = -2*(cross(rob-ro1,l3)).'*cross(n1,l3)-U3p1U3;
kA(2) = -2*(cross(rob-ro2,l3)).'*cross(n2,l3)-U3p2U3;
kA(3) = (-2*(cross(rob-ro3,l3)).'*cross(n3,l3)-U3p3U3);
kA = simplify(kA)
kG = L.'*kA-U.'*D2U*U;
kG = simplify(kG,1000)
%%
val = eig(G);
flag = 0;
if all(val<=1E-9)
    flag = 1;
end


function Ax = CrMat(a)
Ax = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end
function As = sTra(A)
As = 0.5*(A+A.');
end