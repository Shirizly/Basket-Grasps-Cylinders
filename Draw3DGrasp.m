function varargout = Draw3DGrasp(r,h,hcm,t,p1,p2,p3)
fig = figure;
com = [0;0;hcm];
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
n3 = [0,-sin(t),cos(t)].';
ds = abs((p1-p3).'*n3);
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
p1_rel = Rt.'*p1-[0;0;hcm+ds];
phi = atan2(p1_rel(1),p1_rel(2));
n1 = -[sin(phi);cos(phi)*cos(t);cos(phi)*sin(t)];
n2 = n1;
n2(1) = -n2(1);
com = Rt*com;
[X,Y,Z] = cylinder(1,300);
hs = surf(r*X,r*Y,h*Z);
rotate(hs,[1 0 0],rad2deg(t),[0 0 0]);
hs.FaceAlpha = 0.1;
hs.EdgeAlpha = 0.1;
hold on
axis
e1 = [0;0;0];
e2 = [0;-h*sin(t);h*cos(t)];
axisends = [e1 e2];
plot3(axisends(1,:),axisends(2,:),axisends(3,:),'--','lineWidth',3);

pl{1} = plot3(com(1),com(2),com(3),'*r','lineWidth',10,'markerSize',4,'displayName','COM');
% plot3([com(1),com(1)],[com(2),com(2)],[0,com(3)*2],'r--');
pl{2} = plot3(p1(1),p1(2),p1(3),'ok','markerSize',6,'lineWidth',4,'displayName','O_1');
pl{3} = plot3(p2(1),p2(2),p2(3),'ob','markerSize',6,'lineWidth',4,'displayName','O_2');
pl{4} = plot3(p3(1),p3(2),p3(3),'og','markerSize',6,'lineWidth',4,'displayName','O_3');

% plot force lines



f3 = [p3 p3+2*n3];
plot3(f3(1,:),f3(2,:),f3(3,:),'g--','linewidth',3);
f1 = [p1 p1+2*r*n1];
plot3(f1(1,:),f1(2,:),f1(3,:),'k--','linewidth',3);
f2 = [p2 p2+2*r*n2];
plot3(f2(1,:),f2(2,:),f2(3,:),'b--','linewidth',3);
xlabel('X')
ylabel('Y')
zlabel('Z')
legend([pl{:}])
axis equal
% xlim([-1.5,1.5]);
% ylim([-1.8,1]);
% zlim([-0.5,2.5]);

hold off
grid on
% view([-0.5,-0.1,0.4])
switch nargout
    case 1
        varargout{1} = fig;
end
end