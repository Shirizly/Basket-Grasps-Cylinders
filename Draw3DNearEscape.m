function Draw3DNearEscape(r,h,hcm,mincom,t1,t2,p1,p2,p3,col)
hold on
com = [0;0;hcm];
Rx = [1 0 0; 0 cos(t1) -sin(t1); 0 sin(t1) cos(t1)];
Ry = [cos(t2) 0 sin(t2);0 1 0; -sin(t2) 0 cos(t2)];
Rt = Rx*Ry;
com = Rt*com;
[X,Y,Z] = cylinder(1,300);
transl = mincom-com;
hs = surf(r*X+transl(1),r*Y+transl(2),h*Z+transl(3));

rotate(hs,[1 0 0],rad2deg(t1),transl);
rotate(hs,(Rx*[0 1 0].').',rad2deg(t2),transl);
hs.FaceAlpha = 0.01;
hs.EdgeAlpha = 0.05;
hs.FaceColor = col;
hs.EdgeColor = col;

e1 = [0;0;0];
e2 = [0;0;h];
e2 = Rt*e2;
% e2 = [0;-h*sin(t);h*cos(t)];
axisends = [e1+mincom-com e1+mincom-com+e2];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'--');

plot3(mincom(1),mincom(2),mincom(3),'*','Color',col,'markerSize',15,'lineWidth',4)
plot3([mincom(1),mincom(1)],[mincom(2),mincom(2)],[0,mincom(3)*2],'r--');
plot3(p1(1),p1(2),p1(3),'or','markerSize',4,'lineWidth',4);
plot3(p2(1),p2(2),p2(3),'ob','markerSize',4,'lineWidth',4);
plot3(p3(1),p3(2),p3(3),'og','markerSize',4,'lineWidth',4);

% plot force lines


% f1 = [p1 p1+n1];
% plot3(f1(1,:),f1(2,:),f1(3,:),'k--');
% f2 = [p2 p2+n2];
% plot3(f2(1,:),f2(2,:),f2(3,:),'k--');
% f3 = [p3 p3+2*n3];
% plot3(f3(1,:),f3(2,:),f3(3,:),'k--');


xlabel('X')
ylabel('Y')
zlabel('Z')

axis equal
hold off
end