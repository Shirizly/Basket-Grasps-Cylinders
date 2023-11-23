%% Basic setup
syms t ph hcm r ds %ds1
% construct rotation matrix for angle t relative to [0;0;1] in the yz plane
ph1 = ph; ph2 = -ph;
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
% position of the contact points
po1 = [0;0;hcm+ds]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;hcm+ds]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;

g = [0;0;-1];


%% contact normals
n1 = [-sin(ph1),-cos(ph1)*cos(t),-cos(ph1)*sin(t)].';
n2 = [-sin(ph2),-cos(ph2)*cos(t),-cos(ph2)*sin(t)].';
n3 = [0,-sin(t),cos(t)].';
com = [0;0;hcm];
com = Rt*com;
%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = -N\g;

%% finding the generator l as function of ds1,ds2,ph1,ph2
A = [n1(1:2) -n2(1:2)];
dp = p2-p1;
sol = A\dp(1:2);
inter1 = p1+sol(1)*n1;
inter2 = p2+sol(2)*n2;

%% finding p3
h12 = cos(t)*ds+com(3);
h3 = ds*(cos(t)^2+1)/cos(t)+com(3);
s3 = (h3-h12)*sin(t);
p3 = s3*[0;cos(t);sin(t)];

%% setup for basket Condition
l1 = n1; l2 = n2; l3 = n3;
Z1 = inter1; Z2 = inter2; Z3 = [inter1(1:2);h3];
z1 = Z1(3); z2 = Z2(3); z3 = Z3(3);
ro1 = p1; ro2 = p2; ro3 = p3;
psy1 = acos(dot(l1,g)); psy2 = acos(dot(l2,g)); psy3 = acos(dot(l3,g));
zcm = com(3);
c1 = sqrt(1+z1.'*z1*(sin(psy1))^2);
c2 = sqrt(1+z2.'*z2*(sin(psy2))^2);
c3 = sqrt(1+z3.'*z3*(sin(psy3))^2);
p1cond = z1 + L(1)/c1*(z3-z1)*l1;
p2cond = z2 + L(2)/c2*(z3-z2)*l2;
lc1 = L(1);%/c1;
lc2 = L(2);%/c2;
lc3 = L(3);%/c3;

%% condition
u1 = [cross(p2,l1);l1];
u2 = [cross(p1,l2);l2];
u3 = [cross(p1,l3);l3];
star = [com(3) 0 -0.5*com(1);0 com(3) -0.5*com(2); -0.5*com(1) -0.5*com(2) 0];
Sum1 = lc1*sTra(CrMat(ro1).'*CrMat(l1))+lc2*sTra(CrMat(ro2).'*CrMat(l2))+...
    lc3*sTra(CrMat(ro3).'*CrMat(l3));
Gpart1 = [l1 l2 l3].'*(star - Sum1)*[l1 l2 l3];
LB12 = [1/r*eye(2) zeros(2,1); zeros(1,3)];
base_v = [0,1,0].';
v1 = cross(base_v,l1);
Vx1 = CrMat(v1);
cos1 = base_v.'*l1;
s1sq = v1.'*v1;
rot1 = eye(3) + Vx1 + Vx1^2*(1-cos1)/s1sq;
% LB1 = rot1*LB12;
LB1 = LB12;
v2 = cross(base_v,l2);
Vx2 = CrMat(v2);
cos2 = base_v.'*l2;
s2sq = v2.'*v2;
rot2 = eye(3) + Vx2 + Vx2^2*(1-cos2)/s2sq;
% LB2 = rot2*LB12;
LB2 = LB12;
LB3 = zeros(3);

 
Sum21 = lc1*[eye(3) -CrMat(ro1); zeros(3) CrMat(l1)].'*[LB1 -eye(3);-eye(3) zeros(3)]*...
    [eye(3) -CrMat(ro1); zeros(3) CrMat(l1)];
Sum22 = lc2*[eye(3) -CrMat(ro2); zeros(3) CrMat(l2)].'*[LB2 -eye(3);-eye(3) zeros(3)]*...
    [eye(3) -CrMat(ro2); zeros(3) CrMat(l2)];
Sum23= lc3*[eye(3) -CrMat(ro3); zeros(3) CrMat(l3)].'*[LB3 -eye(3);-eye(3) zeros(3)]*...
    [eye(3) -CrMat(ro3); zeros(3) CrMat(l3)];
Sum2 = Sum21 + Sum22 + Sum23;
Gpart2 = [u1 u2 u3].'*Sum2*[u1 u2 u3];
%%
G = simplify(Gpart1 + Gpart2,200)
L = eig(G);


function Ax = CrMat(a)
Ax = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end
function As = sTra(A)
As = 0.5*(A+A.');
end