function [f,h] = draw3DSpheres(r,p1,p2,p3)
%draws a sphere with radius r surrounding each finger
f = figure;
    hold on
    plot3(p1(1),p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
    plot3(p2(1),p2(2),p2(3),'ok','markerSize',4,'lineWidth',4);
    plot3(p3(1),p3(2),p3(3),'ok','markerSize',4,'lineWidth',4);
    [x,y,z] = sphere(30);
    h(1) = surf(r*x+p1(1),r*y+p1(2),r*z+p1(3),'faceAlpha',0.7,'faceColor','y','displayName','S_1,S_2');
    surf(r*x+p2(1),r*y+p2(2),r*z+p2(3),'faceAlpha',0.7,'faceColor','y');
    h(2) = surf(r*x+p3(1),r*y+p3(2),r*z+p3(3),'faceAlpha',0.7,'faceColor','b','displayName','S_3');
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
end