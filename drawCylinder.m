function labeled = drawCylinder(r,h,p,new,col,DN,drawBases,rotation)
if ~new
    figure;
end
[X,Y,Z] = cylinder(r,50);
if ~isempty(rotation)
    translate = p;
    p = zeros(3,1);
end
Z = Z*h+p(3);
Y = Y+p(2);
X = X+p(1);
gr = hgtransform;
% s=surf(X,Y,Z,'parent',gr,'edgeColor','none');
% 
% s.FaceColor=col;
% s.DisplayName = DN;
% 
% s.EdgeColor = col;
% s.EdgeAlpha = 0.1;
hold on
if drawBases
    c = [];
% c(1) = drawCircle(p(1),p(2),p(3),r,col);
c(end+1) = drawCircle(p(1),p(2),p(3)+h,r,col);
set(c,'parent',gr);
end
if ~isempty(rotation)
    
    M = makehgtform('translate',translate,'axisrotate',rotation(1:3),rotation(4));
%     rotate(s,rotation(1:3),rotation(4));
%     if drawBases
%     rotate(c(1),rotation(1:3),rotation(4));
%     M = makehgtform()
    set(gr,'Matrix',M)
    drawnow
%     rotate(c(2),rotation(1:3),rotation(4));
%     end
end
labeled = 0;% s;
end

function circle = drawCircle(xc,yc,zc,r,col)
teta = linspace(0,2*pi,50);
x=xc+r*cos(teta);
y=yc+r*sin(teta) ;
z = zc+zeros(size(x));
circle = patch(x,y,z,col,'edgeColor','none');
end

%% example code for drawing the objects for the paper:
% p=zeros(3,1);
% drawCylinder(0.1,2,p,0)
% p=[6;0;0];
% hold on
% drawCylinder(5,2,p,1)
% axis equal
% p=[12;0;0];
% hold on
% drawCylinder(0.5,2,p,1)
% set(gcf,'renderer','painters')