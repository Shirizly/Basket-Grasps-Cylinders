function Draw3D(r,h,hcm,t,p1,p2,p3,ph1,ph2)
figure;
com = [0;0;hcm];
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
n1 = [-sin(ph1),-cos(ph1)*cos(t),-cos(ph1)*sin(t)].';
n2 = [-sin(ph2),-cos(ph2)*cos(t),-cos(ph2)*sin(t)].';
n3 = [0,-sin(t),cos(t)].';

com = Rt*com;
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
plot3(p1(1),p1(2),p1(3),'or','markerSize',4,'lineWidth',4);
plot3(p2(1),p2(2),p2(3),'ob','markerSize',4,'lineWidth',4);
plot3(p3(1),p3(2),p3(3),'og','markerSize',4,'lineWidth',4);

% plot force lines


f1 = [p1 p1+n1];
plot3(f1(1,:),f1(2,:),f1(3,:),'k--');
f2 = [p2 p2+n2];
plot3(f2(1,:),f2(2,:),f2(3,:),'k--');
f3 = [p3 p3+2*n3];
plot3(f3(1,:),f3(2,:),f3(3,:),'k--');


xlabel('X')
ylabel('Y')
zlabel('Z')

axis equal
hold off
end