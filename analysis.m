% syms t p1 p2
t = pi()/5;
ph = pi()*2/3;
ph1 = ph; ph2 = -ph;
r = 0.5;
h = 2;
hcm = 0.5*h;
com = [0;0;hcm];
ds = 0.4;
ds1 = ds;
ds2 = ds;
s1 = r + hcm + ds1;
s2 = r + hcm + ds2;
% construct rotation matrix for angle t relative to [0;0;1] in the yz plane
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
% position of the contact points
po1 = [0;0;s1-r]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;s2-r]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;
% contact normals
n1 = [-sin(ph1),-cos(ph1)*cos(t),-cos(ph1)*sin(t)].';
n2 = [-sin(ph2),-cos(ph2)*cos(t),-cos(ph2)*sin(t)].';
n3 = [0,-sin(t),cos(t)].';

com = Rt*com;


%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = -N\g;

%% finding p3 for coplanar p1p2
%  syms ds h3
if ds1 == ds2 &&ph1==-ph2
    h12 = cos(t)*ds;
    h3 = ds*(cos(t)^2+1)/cos(t);
    s3 = (h3-h12)*sin(t);
    if s3>0 && s3<r
        p3 = s3*[0;cos(t);sin(t)];
    else
        p3 = [];
    end
else
    A = [n1(1:2) -n2(1:2)];
    dp = p2-p1;
    sol = A\dp(1:2);
    inter1 = p1+sol(1)*n1;
    % inter2 = p2+sol(2)*n2;
    %     h1 = inter1(3)-com(3);
    %     h2 = inter2(3)-com(3);
    lxy = inter1(1:2)-com(1:2);
    h1 = -(ds1*cos(ph1)*sin(ph2) - ds2*cos(ph1)*sin(ph2) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
    h2 = -(ds1*cos(ph2)*sin(ph1) - ds2*cos(ph2)*sin(ph1) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
    h3 = -(ds1*cos(ph1)*sin(ph2) - ds2*cos(ph2)*sin(ph1) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
    inter3 = [lxy;h3]+com;
    dist3 = n3.'*inter3;
    p3 = inter3-n3*dist3;
end
 
%% plotting the contacts and force lines
figure(1);
subplot(1,2,1) % yz plane
hold on

%plot rectangle view of the cylinder
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
plot(com(2),com(3),'*r');
plot([com(2),com(2)],[0,com(3)*2],'r--');
%plot contact points
plot(p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
plot(p2(2),p2(3),'ok','markerSize',4,'lineWidth',4);
plot(p3(2),p3(3),'ok','markerSize',4,'lineWidth',4);
% plot force lines
f1 = [p1 p1+n1];
plot(f1(2,:),f1(3,:),'k--');
f2 = [p2 p2+n2];
plot(f2(2,:),f2(3,:),'k--');
f3 = [p3 p3+2*n3];
plot(f3(2,:),f3(3,:),'k--');
xlabel('Y')
ylabel('Z')

axis equal
hold off

%
subplot(1,2,2) % yx plane
hold on

% plot axis of cylinder
e1 = [0;0];
e2 = [-h*sin(t);0];
axisends = [e1 e2];
plot(axisends(1,:),axisends(2,:),'b--');

% mark com
plot(com(2),com(1),'*r')

%plot contact points
plot(p1(2),p1(1),'ok','markerSize',4,'lineWidth',4)
plot(p2(2),p2(1),'ok','markerSize',4,'lineWidth',4)
plot(p3(2),p3(1),'ok','markerSize',4,'lineWidth',4)

% plot force lines
f1 = [p1 p1+n1];
plot(f1(2,:),f1(1,:),'k--');
f2 = [p2 p2+n2];
plot(f2(2,:),f2(1,:),'k--');
f3 = [p3 p3+2*n3];
plot(f3(2,:),f3(1,:),'k--');
xlabel('Y')
ylabel('X')
hold off
axis equal
%% checking stability
[flag, eigv] = MasonConditionFunction(r,hcm,t,ph,ds)
[flag,val] = MinCondition2D(hcm,r,ph,t,ds)

%% plot in 3-D
figure(2);
[X,Y,Z] = cylinder(1,300);
hs = surf(r*X,r*Y,h*Z);
rotate(hs,[1 0 0],rad2deg(t),[0 0 0]);
hs.FaceAlpha = 0.01;
hs.EdgeAlpha = 0.05;
hold on
axis
e1 = [0;0;0];
e2 = [0;-h*sin(t);h*cos(t)];
axisends = [e1 e2];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'--');

plot3(com(1),com(2),com(3),'*r')
plot3([com(1),com(1)],[com(2),com(2)],[0,com(3)*2],'r--');
plot3(p1(1),p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
plot3(p2(1),p2(2),p2(3),'ok','markerSize',4,'lineWidth',4);
plot3(p3(1),p3(2),p3(3),'ok','markerSize',4,'lineWidth',4);

% plot force lines
f1 = [p1 p1+n1];
plot3(f1(1,:),f1(2,:),f1(3,:),'k--');
f2 = [p2 p2+n2];
plot3(f2(1,:),f2(2,:),f2(3,:),'k--');
f3 = [p3 p3+2*n3];
plot3(f3(1,:),f3(2,:),f3(3,:),'k--');

% plot generators (projection to 3-D) (for coplanar p1p2)
g1p = (s1-r)*[0;-sin(t);cos(t)];
g1d = [0;0;-h/2];
g1 = [g1p-g1d g1p+g1d];
plot3(g1(1,:),g1(2,:),g1(3,:),'g','lineWidth',4);
g2p = (s1-r)*[0;-sin(t);cos(t)];
g2d = [0;r;0];
g2 = [g2p-g2d g2p+g2d];
plot3(g2(1,:),g2(2,:),g2(3,:),'g','lineWidth',4);
g3p = (s1-r)*[0;-sin(t);cos(t)]+s3*[0;cos(t);sin(t)];
g3d = [r;0;0];
g3 = [g3p-g3d g3p+g3d];
plot3(g3(1,:),g3(2,:),g3(3,:),'g','lineWidth',4);
xlabel('X')
ylabel('Y')
zlabel('Z')

axis equal
hold off