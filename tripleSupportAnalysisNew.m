syms r t1 ds phi hcm
ph1 = phi; ph2 = -phi;
vx = [1;0;0];
vy = [0;cos(t1);sin(t1)];
vz = [0;-sin(t1);cos(t1)];
Rt = [vx vy vz];
com_0 = [0;0;hcm];
com = Rt*com_0;

% position of the contact points
po1 = [0;0;hcm+ds]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;hcm+ds]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;

p3 = ds*tan(t1)*[0;cos(t1);sin(t1)];
%%
syms r p1x p1y p1z p3y p3z
p1 = [p1x;p1y;p1z];
p2 = p1;
p2(1) = -p2(1);
p3 = [0;p3y;p3z];

x = sym('x',[2,1]);
t1t = x(1);
t2t = x(2);

p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
syms phi3 d3
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rt*Rz*[0 -d3 0].';
n = Rt*nt;
p1p = p1t-(p1t.'*n)*n;
p2p = p2t-(p2t.'*n)*n;
% rad1 = p1p-center;
rad1 = cross(p1t-center,n);
% rad2 = p2p-center;
rad2 = cross(p2t-center,n);
eq1 = rad1.'*rad1-r^2;
eq2 = rad2.'*rad2-r^2;
eq = [eq1;eq2];
eq = simplify(eq)
pretty(eq)
%%
eqex = expand(eq);
eq3 = simplify(eq(1)-eq(2))
% t1_sol = solve(eq3,x(1))
t2_sol_withCond = solve(eq3,x(2), 'Real', true)
t2_sol = simplify(t2_sol_withCond)
% t2_sol_cond = t2_sol_withCond.conditions
%%
% eq1 = simplify(subs(eq(1),x(2),t2_sol(2)));
% eq2 = simplify(subs(eq(1),x(2),t2_sol(3)));
eq1 = (cos(x(1))*(p1x - d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) + (d3*sin(phi3)*(p1z - p3z + d3*cos(phi3)*sin(x(1)) - (d3^2*cos(x(1))*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2 - r^2 - ((d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2 - 1)*(d3*cos(phi3) + p1y*cos(x(1)) - p3y*cos(x(1)) + p1z*sin(x(1)) - p3z*sin(x(1)))^2 + (sin(x(1))*(p1x - d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) - (d3*sin(phi3)*(p1y - p3y + d3*cos(phi3)*cos(x(1)) + (d3^2*sin(phi3)^2*sin(x(1)))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2;
 
eq2 = (cos(x(1))*(p1x + d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) - (d3*sin(phi3)*(p1z - p3z + d3*cos(phi3)*sin(x(1)) - (d3^2*cos(x(1))*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2 - r^2 - ((d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2 - 1)*(d3*cos(phi3) + p1y*cos(x(1)) - p3y*cos(x(1)) + p1z*sin(x(1)) - p3z*sin(x(1)))^2 + (sin(x(1))*(p1x + d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) + (d3*sin(phi3)*(p1y - p3y + d3*cos(phi3)*cos(x(1)) + (d3^2*sin(phi3)^2*sin(x(1)))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2;
pretty(simplify(eq1))
pretty(simplify(eq2))

% eq4_1 = simplify(subs(eq(1),x(2),t2_sol(1)))
% pretty(eq4_1)
% t1_sol1 = simplify(solve(eq4_1,x(1)))
% %%
% eq4_2 = simplify(subs(eq(1),[x(2)],[t2_sol(2)]));
% syms A psi
% eq4_2 = simplify(subs(eq4_2,p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)),A*sin(x(1)+psi)))
% pretty(eq4_2)
% %%
% t1_sol2 = simplify(solve(eq4_2,x(1)), 'Real', true)
% %%
% eq4_3 = simplify(subs(eq(1),[x(2),k],[t2_sol(3),0]))
% pretty(eq4_3)
% t1_sol3 = simplify(solve(eq4_3,x(1)), 'Real', true)