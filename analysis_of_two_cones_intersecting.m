
p1n = sym('p1',[3,1]);
p2n = sym('p2',[3,1]);

syms r

syms x y z ph1 ph2

cone = [x;y;z];

eq1 = cone.'*p1n-cos(ph1);
eq2 = cone.'*p2n-cos(ph2);
eq3 = x^2+y^2+z^2-1

sol = solve([eq1;eq2;eq3],[x;y;z]);
solv = [sol.x sol.y sol.z];
% matlabFunction(solv,'file','comGivenBaseCenter')% don't overwrite for no
% reason, made changes in file to simplify

