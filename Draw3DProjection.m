function fig = Draw3DProjection(p1,p2,p3,h,r,t,p1t,p2t)
fig = figure();
[X,Y,Z] = cylinder(1,300);
X = X-p3(1);
Y = Y-p3(1);
Z = Z-p3(3);
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

plot3(p1t(1),p1t(2),p1t(3),'or','markerSize',6,'lineWidth',4);
plot3(p2t(1),p2t(2),p2t(3),'or','markerSize',6,'lineWidth',4);
plot3(p3t(1),p3t(2),p3t(3),'or','markerSize',6,'lineWidth',4);
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