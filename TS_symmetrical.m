function [t1s1,t1s2] = TS_symmetrical(r,p1,p2,p3,phi)
% Function returns the theta1 values at the object stances that maintain
% triple contact, with the finger base at the base edge, and the stance
% symmetric relative to the YZ-plane

% This computation is based on a projection of the fingers O1,O2 to the
% YZ-plane
% start by finding p12 which represent p1,p2 in the projected plane.
p12 = (p1+p2)/2;

% %extract phi1=phi2 from finger positions and r - not necessary anymore
% adp = norm(p1-p2)/2; %half distance between p1,p2
% bdp = sqrt(r^2-adp^2); %other edge of the right triangle
% phi = pi()-atan2(adp,bdp);
%% first case, p3 touches farther 'vertex' of the rectangle from the edge supported by p12
% compute angle between line connecting p3 and average p12, and the base of
% the cylinder, projected to y-z plane, use to compute t1 for first case
d1 = r*(1-cos(phi));
dp = p12-p3; %distance between p12 and p3
dpnorm = norm(dp);
if d1<=dpnorm %distance between fingers allow for this stance to exist
    t1s1 = asin(d1/dpnorm)+atan2(dp(2),dp(3));
else %fingers too close together
    t1s1 = [];
end
%% second case, p3 touches 'vertex' of the edge supported by p12
d1 = r*(1+cos(phi));
t1s2 = -(asin(d1/dpnorm)-atan2(dp(2),dp(3)));
% note that it is possible for this stance to require a cylinder taller
% than the actual object is, but that is checked outside this function
end