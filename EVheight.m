function sol = EVheight(dp,r,hcm)
% receive vector between fingers, cylinder parameters. check first for each
% possibile stance if fingers are not too far apart (or too close for case
% 1)
sol1 = inf*ones(3,1);
sol2 = inf*ones(3,1);
% case 1, base finger on farther vertex of rectangle.
if sqrt(dp.'*dp-4*r^2)<=2*hcm && norm(dp)>=2*r %checking that the fingers are not too far apart
    theta1_1 = pi()/2-atan2(abs(dp(2)),abs(dp(1)))-asin(2*r/norm(dp));
    sol1 = r*[cos(theta1_1);sin(theta1_1)] + hcm*[sin(theta1_1);cos(theta1_1)];
end
% case 2, base finger on closer vertex of rectangle (essentially two
% fingers on same edge).
if norm(dp)<2*hcm %checking that the fingers are not too far apart
theta1_2 = asin(abs(dp(1)/norm(dp)));
sol2 = hcm*[sin(theta1_2);cos(theta1_2)];
end
if all([sol1;sol2]==inf)
    sol = [sol1;NaN];
    return
end
if sol1(2)<sol2(2) % return the stance with lower com.
    sol = [sol1;theta1_1];
else
    sol = [sol2;theta1_2];
end
end