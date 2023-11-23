f = figure;
    hold on
    plot3(p1(1),p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
    plot3(p2(1),p2(2),p2(3),'ob','markerSize',4,'lineWidth',4);
    plot3(p3(1),p3(2),p3(3),'og','markerSize',4,'lineWidth',4);
    [x,y,z] = sphere;
    mesh(r*x+p1(1),r*y+p1(2),r*z+p1(3),'faceAlpha',0.2);
    mesh(r*x+p2(1),r*y+p2(2),r*z+p2(3),'faceAlpha',0.2);
%     mesh(r*x+p3(1),r*y+p3(2),r*z+p3(3),'faceAlpha',0.2);
c = centersol(20,:).';
[X,Y,Z]=cylinder([0 .5],50 );
% axis([0 1,-1 1,-.5 .5])
% M=makehgtform('translate',[0,0,0],'xrotate',0,'yrotate',0);
hs=surf(X+c(1),Y+c(2),Z+c(3),'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.4);

axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')