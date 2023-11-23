function sol = doubleSupportEscape(r,hcm,p1,p2,p3)
% function receives cylinder parameters r and hcm, and physical position of
% fingers. Compute all double-support escape stances and return the minimum
% of their height.
% for now I build the function to analyze only the isosceles finger
% formation
solrectangle = rectangleDSHeight(r,hcm,p1,p3);
heightRectangle = solrectangle(2);
heightCircle = circleDSHeight(r,p1,p2);
sol = [heightRectangle,heightCircle];
end

function sol = rectangleDSHeight_old(r,hcm,p1,p3)
[f1,f3] = RZFromXYZ(p1,p3);
% first - find the orientation of the rectangle in edge-edge stance
% (using solution of 2-D EE-escape)
% requires defining the problem as EE support on two edges
dp = f1-f3; %inter-finger vector
nj = [0;1]; %normals for the edges
nk = [1;0];
vj = [-r;-hcm];
vk = [-r;hcm];
dels1 = 2*r;
dels2 = 2*hcm;
sol = EEheight(dp,nj,nk,vj,vk,dels1,dels2);
% receive multiple solutions and need to identify the relevant one
if size(sol,2)>1
    disp('multiple edge-edge solutions');
end
if size(sol,2)==1
    sol(1:2) = sol(1:2) + f3;
    return
end
% if no solution was found, that means the fingers can't support
% the object in a two finger edge-edge saddle, and instead need to consider
% edge-vertex. It is unclear at this point if it is possible for
% double-support edge-vertex stances to be lower than the lowest
% triple-support stance, requires checking.
sol = EVheight(dp,r,hcm);
sol(1:2) = sol(1:2) + f3;
end

function sol = rectangleDSHeight(r,hcm,p1,p3)
[p_c1,p_c3] = RZFromXYZ(p1,p3);
% first - find the orientation of the rectangle in edge-edge stance
% (using solution of 2-D EE-escape)
% requires defining the problem as EE support on two edges

% vj = [r;-hcm];
% vk = [-r;-hcm];
dels1 = 2*r;
dels2 = 2*hcm;
sol = EEheight_new(r,hcm,p_c1,p_c3,dels1,dels2);

if size(sol,2)==1
    sol(1:2) = sol(1:2) + p_c3;
    return
end
sol = [0;inf];
end

function height2 = circleDSHeight(r,p1,p2)
[f1,f2] = RZFromXYZ(p1,p2);
% the height is taken relative to the fingers, using Pyhagoras given a
% finger and the center of the circle as two endpoints of the hypoten.
heightRelToFingers = sqrt(r^2-(norm(f2-f1)/2)^2);
height2 = f1(2) + heightRelToFingers;
end

function [f1,f2] = RZFromXYZ(p1,p2)
%takes an XYZ 3-D vector, project both on the plane holding both points and
%gravity
f1 = p1([1,3]);
f12= p2-p1;
f2 = f1;
f2(1) = f1(1) + norm(f12(1:2));
f2(2) = p2(3);
end

