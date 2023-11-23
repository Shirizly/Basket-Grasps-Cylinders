syms r t d psy b phi
% p1 = sym('p1',[3,1]);
% p2 = p1;
% p2(1) = -p2(1);
% p3 = 
hcm = 1;
ph1 = phi;
ph2 = -phi;
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
po1 = [0;0;hcm+b]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;hcm+ds]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;
p3 = b*tan(t)*[0;cos(t);sin(t)];
noon = [0;cos(t);sin(t)];
perpToBoth = [-1;0;0];
clockArm = d*(noon*cos(psy) + perpToBoth*sin(psy));
centerb = clockArm;

d1 = p1-centerb;
nd1 = (d1.'*d1)^0.5;
d1n = d1/nd1;
phi1 = asin(r/nd1);

d2 = p2-centerb;
nd2 = (d2.'*d2)^0.5;
d2n = d2/nd2;
phi2 = asin(r/nd2);

syms x y z

cone = [x;y;z];

eq1 = simplify(cone.'*d1n-cos(phi1));
eq2 = simplify(cone.'*d2n-cos(phi2));
eq3 = x^2+y^2+z^2-1;

sol = solve([eq1;eq2;eq3],[x;y;z]);
solv = [sol.x sol.y sol.z];

com = simplify(solv(2,:).');
%%
corrected_centerb = -(centerb-p3).'*com/(com.'*com)*com+centerb;
corrected_com = com+corrected_centerb;


%%

d2U = hessian(corrected_com(3),[d,psy]);
%%
d2UGrasp = subs(d2U,[d,psy],[0,0]);
%%
matlabFunction(d2UGrasp,'file','HessianOfGrasp')