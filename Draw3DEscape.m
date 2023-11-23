function Draw3DEscape(r,h,hcm,mincom,t1,t2,phi3,p1,p2,p3)
hold on
com = [0;0;hcm];
Rx = [1 0 0; 0 cos(t1) -sin(t1); 0 sin(t1) cos(t1)];
Ry = [cos(t2) 0 sin(t2);0 1 0; -sin(t2) 0 cos(t2)];
Rt = Rx*Ry;
% Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
Rz = eye(3);
com = Rt*com;
[X,Y,Z] = cylinder(1,300);
transl = mincom-com;
hs = surf(r*X+transl(1),r*Y+transl(2),h*Z+transl(3));
% rotate(hs,[0 0 1],rad2deg(phi3),transl);
rotate(hs,[1 0 0],rad2deg(t1),transl);
rotate(hs,(Rx*[0 1 0].').',rad2deg(t2),transl);

hs.FaceAlpha = 0.1;
hs.EdgeAlpha = 0.1;


e1 = [0;0;0];
e2 = [0;0;h];
e2 = Rt*e2;
% e2 = [0;-h*sin(t);h*cos(t)];
axisends = [e1+mincom-com e1+mincom-com+e2];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'--','lineWidth',3);

pl{1} = plot3(mincom(1),mincom(2),mincom(3),'*r','lineWidth',10,'markerSize',4,'displayName','COM');
% plot3([com(1),com(1)],[com(2),com(2)],[0,com(3)*2],'r--');
pl{2} = plot3(p1(1),p1(2),p1(3),'ok','markerSize',6,'lineWidth',4,'displayName','O_1');
pl{3} = plot3(p2(1),p2(2),p2(3),'ob','markerSize',6,'lineWidth',4,'displayName','O_2');
pl{4} = plot3(p3(1),p3(2),p3(3),'og','markerSize',6,'lineWidth',4,'displayName','O_3');

% plot force lines
if norm(p3-transl)<=r
[N,L,~] = TSNormals(mincom,p1,p2,p3,[t1,t2]);
n1 = N(:,1);
n2 = N(:,2);
n3 = N(:,3);
f3 = [p3 p3+2*n3];
plot3(f3(1,:),f3(2,:),f3(3,:),'g--','linewidth',3);
f1 = [p1 p1+2*r*n1];
plot3(f1(1,:),f1(2,:),f1(3,:),'k--','linewidth',3);
f2 = [p2 p2+2*r*n2];
plot3(f2(1,:),f2(2,:),f2(3,:),'b--','linewidth',3);
% f1 = [p1 p1+2*r*n1];
% plot3(f1(1,:),f1(2,:),f1(3,:),'k--');
% f2 = [p2 p2+2*r*n2];
% plot3(f2(1,:),f2(2,:),f2(3,:),'k--');
% f3 = [p3 p3+2*n3];
% plot3(f3(1,:),f3(2,:),f3(3,:),'k--');
tor1 = cross(p1,n1*L(1));
tor2 = cross(p2,n2*L(2));
tor3 = cross(p3,n3*L(3));
torg = cross(mincom,[0;0;-1]);
sumtorque = tor1+tor2+tor3+torg
sumforce = n1*L(1)+n2*L(2)+n3*L(3)
forcecond = all(L>0)
end

xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
% xlim([-1.5,1.5]);
% ylim([-1.8,1]);
% zlim([-0.5,2.5]);
legend([pl{:}])
hold off
grid on
% view([-0.5,-0.1,0.4])
end