syms t ph1  hcm r base height ph2 ds2 %ds1
% construct rotation matrix for angle t relative to [0;0;1] in the yz plane
ds1 = 0;
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
% position of the contact points
po1 = [0;0;hcm+ds1]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;hcm+ds2]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;
h_p12 = simplify(p1(3))
eq1_horizontal = p2(3)==p1(3);
ds2 = solve(eq1_horizontal,ds2)
p2 = simplify(subs(p2))
% bias = p1(3)*tan(t);
% a = r/cos(t); %ellipse width/2
% b = r; %ellipse height/2
% eq1_ellipse = p2(1)^2/a^2+(p2(2)-bias)^2/b^2 == 1; 
eq2_baselength = (p2(1:2)-p1(1:2)).'*(p2(1:2)-p1(1:2)) == base^2
sol = solve(eq2_baselength,ph2,'ReturnConditions',true)

%% find p3 based on parametric p1 and p2
syms p1x p1y p1z p2x p2y p2z p3x p3y p3z height slope
p1 = [p1x;p1y;p1z];
p2 = [p2x;p2y;p2z];
p3 = [p3x;p3y;p3z];
base_vector = p1-p2;
base_middle = (p1+p2)/2;
v2p3 = (p3-base_middle);
eq1 = v2p3.'*base_vector==0;
eq2 = v2p3.'*v2p3==height^2;
eq3 = v2p3(3)/(v2p3(1:2).'*v2p3(1:2))^0.5==slope;
p3v = solve([eq1,eq2,eq3],p3)
p3v1 = [p3v.p3x(1);p3v.p3y(1);p3v.p3z(1)]
p3v2 = [p3v.p3x(2);p3v.p3y(2);p3v.p3z(2)]
p3v3 = [p3v.p3x(3);p3v.p3y(3);p3v.p3z(3)]
p3v4 = [p3v.p3x(4);p3v.p3y(4);p3v.p3z(4)]

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
%  syms ds h3
% if ds1 == ds2 % for coplanar p1p2 (3.a)
%     h12 = cos(t)*ds;
%     h3 = ds*(cos(t)^2+1)/cos(t);
%     s3 = (h3-h12)*sin(t);
%     if s3>0 && s3<r
%         p3 = s3*[0;cos(t);sin(t)];
%     else
%         p3 = [];
%     end
% else % for general case (3.b)
    h1 = inter1(3)-com(3);
    h2 = inter2(3)-com(3);
    lxy = simplify(inter1(1:2)-com(1:2));
    h3y = (lxy(2)-h1*L(1)*n1(2)-h2*L(2)*n2(2))/(L(3)*n3(2));
% end
h3 = -(ds1*cos(ph1)*sin(ph2) - ds2*cos(ph2)*sin(ph1) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
inter3 = [lxy;h3]+com;
dist3 = n3.'*inter3;
p3 = inter3-n3*dist3;
syms base height slope
middle_base = (p1+p2)/2;
eq1 = norm(middle_base-p3)==height;
triangle_normal = middle_base-p3;
eq2 = (p2-p1).'*triangle_normal==0;
sol = solve([eq1,eq2],[ph2,ds2])
