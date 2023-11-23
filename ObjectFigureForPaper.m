p=[3;0;2];
rotation = [0;1;0;pi()/2];
labeled(1) = drawCylinder(0.35,7,p,0,'r','example 1',1,rotation);

rotation = [];
p=[5;0;0];
hold on
labeled(2) = drawCylinder(2,0.8,p,1,'g','example 2',1,rotation);
p=[12;0;0];
labeled(3) = drawCylinder(1,4,p,1,'b','example 3',1,rotation);
axis equal
set(gcf,'renderer','painters')
light
set(gca,'visible','off')
view([2,-10,2])
% legend(labeled)

%%
epsclean('caps.eps','capsClean.eps')