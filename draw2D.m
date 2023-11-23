function draw2D(r,t,h,com,p1,p2,p3,n1,n2,n3)
figure();
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
end