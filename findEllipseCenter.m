function c = findEllipseCenter(a,b,th,p1,p2)
pc = (p1+p2)/2;
R = [cos(th) sin(th);-sin(th) cos(th)];
p1 = pc+R*(p1-pc);
p2 = pc+R*(p2-pc);
options = optimoptions('fsolve','Display','none');
f = @(c)[1/a^2*(p1(1)^2-2*p1(1)*c(1)+c(1)^2)+1/b^2*(p1(2)^2-2*p1(2)*c(2)+c(2)^2)-1;...
1/a^2*(p2(1)^2-2*p2(1)*c(1)+c(1)^2)+1/b^2*(p2(2)^2-2*p2(2)*c(2)+c(2)^2)-1];

c = fsolve(f,[0;0],options);
c = pc+R.'*(c-pc);
end

function c = findEllipseCenterNoAngle(a,b,p1,p2)
options = optimoptions('fsolve','Display','none');
f = @(c)[1/a^2*(p1(1)^2-2*p1(1)*c(1)+c(1)^2)+1/b^2*(p1(2)^2-2*p1(2)*c(2)+c(2)^2)-1;...
1/a^2*(p2(1)^2-2*p2(1)*c(1)+c(1)^2)+1/b^2*(p2(2)^2-2*p2(2)*c(2)+c(2)^2)-1];

c = fsolve(f,[0;0],options);
end
%% didn't work - norm "vanishes" effect of rotation
% function c = findEllipseCenter(a,b,th,p1,p2)
% options = optimoptions('fsolve','Display','none');
% f = @(c)[((cos(th)*(c(1) - p1(1)))/a - (sin(th)*(c(2) - p1(2)))/b)^2 + ((cos(th)*(c(2) - p1(2)))/b + (sin(th)*(c(1) - p1(1)))/a)^2-1;
% ((cos(th)*(c(1) - p2(1)))/a - (sin(th)*(c(2) - p2(2)))/b)^2 + ((cos(th)*(c(2) - p2(2)))/b + (sin(th)*(c(1) - p2(1)))/a)^2-1];
% 
% c = fsolve(f,[0;0],options);
% end