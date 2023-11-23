syms z ph1 t base r h
% r^2*z1^2 + r^2*exp(ph1*2i) + r^2*exp(ph1*2i + t*4i) + 4*r^2*z1*exp(ph1*1i + t*2i) + 12*r^2*z1*exp(ph1*3i + t*2i) + 2*r^2*z1*exp(ph1*3i + t*4i) + 4*base^2*z1^2*exp(ph1*2i) + 2*r^2*z1^3*exp(ph1*1i) + r^2*z1^2*exp(ph1*4i) + r^2*z1^4*exp(ph1*2i) + r^2*z1^2*exp(t*4i) + 8*base^2*z1^2*exp(ph1*2i + t*2i) + 4*base^2*z1^2*exp(ph1*2i + t*4i) + 12*r^2*z1^3*exp(ph1*1i + t*2i) + 2*r^2*z1^3*exp(ph1*1i + t*4i) + 4*r^2*z1^3*exp(ph1*3i + t*2i) + r^2*z1^2*exp(ph1*4i + t*4i) + r^2*z1^4*exp(ph1*2i + t*4i) + 2*r^2*z1*exp(ph1*3i) ~= 2*r^2*exp(ph1*2i + t*2i) + 2*r^2*z1*exp(ph1*1i + t*4i) + 4*r^2*z1^2*exp(ph1*2i) + 2*r^2*z1^3*exp(ph1*3i) + 2*r^2*z1^2*exp(t*2i) + 24*r^2*z1^2*exp(ph1*2i + t*2i) + 4*r^2*z1^2*exp(ph1*2i + t*4i) + 2*r^2*z1^2*exp(ph1*4i + t*2i) + 2*r^2*z1^4*exp(ph1*2i + t*2i) + 2*r^2*z1^3*exp(ph1*3i + t*4i) + 2*r^2*z1*exp(ph1*1i)
z1 = [root(2*r^2*z^4*exp(ph1*2i)*exp(t*2i) - r^2*z^4*exp(ph1*2i)*exp(t*4i) - r^2*z^4*exp(ph1*2i) - 12*r^2*z^3*exp(ph1*1i)*exp(t*2i) + 2*r^2*z^3*exp(ph1*3i)*exp(t*4i) - 4*r^2*z^3*exp(ph1*3i)*exp(t*2i) - 2*r^2*z^3*exp(ph1*1i)*exp(t*4i) + 2*r^2*z^3*exp(ph1*3i) - 2*r^2*z^3*exp(ph1*1i) + 24*r^2*z^2*exp(ph1*2i)*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i)*exp(t*4i) + 2*r^2*z^2*exp(ph1*4i)*exp(t*2i) - 8*base^2*z^2*exp(ph1*2i)*exp(t*2i) - 4*base^2*z^2*exp(ph1*2i)*exp(t*4i) - r^2*z^2*exp(ph1*4i)*exp(t*4i) + 2*r^2*z^2*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i) - 4*base^2*z^2*exp(ph1*2i) - r^2*z^2*exp(t*4i) - r^2*z^2*exp(ph1*4i) - r^2*z^2 - 12*r^2*z*exp(ph1*3i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i)*exp(t*4i) + 2*r^2*z*exp(ph1*1i)*exp(t*4i) - 4*r^2*z*exp(ph1*1i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i) + 2*r^2*z*exp(ph1*1i) + 2*r^2*exp(ph1*2i)*exp(t*2i) - r^2*exp(ph1*2i)*exp(t*4i) - r^2*exp(ph1*2i), z, 1);...
 root(2*r^2*z^4*exp(ph1*2i)*exp(t*2i) - r^2*z^4*exp(ph1*2i)*exp(t*4i) - r^2*z^4*exp(ph1*2i) - 12*r^2*z^3*exp(ph1*1i)*exp(t*2i) + 2*r^2*z^3*exp(ph1*3i)*exp(t*4i) - 4*r^2*z^3*exp(ph1*3i)*exp(t*2i) - 2*r^2*z^3*exp(ph1*1i)*exp(t*4i) + 2*r^2*z^3*exp(ph1*3i) - 2*r^2*z^3*exp(ph1*1i) + 24*r^2*z^2*exp(ph1*2i)*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i)*exp(t*4i) + 2*r^2*z^2*exp(ph1*4i)*exp(t*2i) - 8*base^2*z^2*exp(ph1*2i)*exp(t*2i) - 4*base^2*z^2*exp(ph1*2i)*exp(t*4i) - r^2*z^2*exp(ph1*4i)*exp(t*4i) + 2*r^2*z^2*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i) - 4*base^2*z^2*exp(ph1*2i) - r^2*z^2*exp(t*4i) - r^2*z^2*exp(ph1*4i) - r^2*z^2 - 12*r^2*z*exp(ph1*3i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i)*exp(t*4i) + 2*r^2*z*exp(ph1*1i)*exp(t*4i) - 4*r^2*z*exp(ph1*1i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i) + 2*r^2*z*exp(ph1*1i) + 2*r^2*exp(ph1*2i)*exp(t*2i) - r^2*exp(ph1*2i)*exp(t*4i) - r^2*exp(ph1*2i), z, 2);...
 root(2*r^2*z^4*exp(ph1*2i)*exp(t*2i) - r^2*z^4*exp(ph1*2i)*exp(t*4i) - r^2*z^4*exp(ph1*2i) - 12*r^2*z^3*exp(ph1*1i)*exp(t*2i) + 2*r^2*z^3*exp(ph1*3i)*exp(t*4i) - 4*r^2*z^3*exp(ph1*3i)*exp(t*2i) - 2*r^2*z^3*exp(ph1*1i)*exp(t*4i) + 2*r^2*z^3*exp(ph1*3i) - 2*r^2*z^3*exp(ph1*1i) + 24*r^2*z^2*exp(ph1*2i)*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i)*exp(t*4i) + 2*r^2*z^2*exp(ph1*4i)*exp(t*2i) - 8*base^2*z^2*exp(ph1*2i)*exp(t*2i) - 4*base^2*z^2*exp(ph1*2i)*exp(t*4i) - r^2*z^2*exp(ph1*4i)*exp(t*4i) + 2*r^2*z^2*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i) - 4*base^2*z^2*exp(ph1*2i) - r^2*z^2*exp(t*4i) - r^2*z^2*exp(ph1*4i) - r^2*z^2 - 12*r^2*z*exp(ph1*3i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i)*exp(t*4i) + 2*r^2*z*exp(ph1*1i)*exp(t*4i) - 4*r^2*z*exp(ph1*1i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i) + 2*r^2*z*exp(ph1*1i) + 2*r^2*exp(ph1*2i)*exp(t*2i) - r^2*exp(ph1*2i)*exp(t*4i) - r^2*exp(ph1*2i), z, 3);...
 root(2*r^2*z^4*exp(ph1*2i)*exp(t*2i) - r^2*z^4*exp(ph1*2i)*exp(t*4i) - r^2*z^4*exp(ph1*2i) - 12*r^2*z^3*exp(ph1*1i)*exp(t*2i) + 2*r^2*z^3*exp(ph1*3i)*exp(t*4i) - 4*r^2*z^3*exp(ph1*3i)*exp(t*2i) - 2*r^2*z^3*exp(ph1*1i)*exp(t*4i) + 2*r^2*z^3*exp(ph1*3i) - 2*r^2*z^3*exp(ph1*1i) + 24*r^2*z^2*exp(ph1*2i)*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i)*exp(t*4i) + 2*r^2*z^2*exp(ph1*4i)*exp(t*2i) - 8*base^2*z^2*exp(ph1*2i)*exp(t*2i) - 4*base^2*z^2*exp(ph1*2i)*exp(t*4i) - r^2*z^2*exp(ph1*4i)*exp(t*4i) + 2*r^2*z^2*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i) - 4*base^2*z^2*exp(ph1*2i) - r^2*z^2*exp(t*4i) - r^2*z^2*exp(ph1*4i) - r^2*z^2 - 12*r^2*z*exp(ph1*3i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i)*exp(t*4i) + 2*r^2*z*exp(ph1*1i)*exp(t*4i) - 4*r^2*z*exp(ph1*1i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i) + 2*r^2*z*exp(ph1*1i) + 2*r^2*exp(ph1*2i)*exp(t*2i) - r^2*exp(ph1*2i)*exp(t*4i) - r^2*exp(ph1*2i), z, 4)];
eq = 2*r^2*z^4*exp(ph1*2i)*exp(t*2i) - r^2*z^4*exp(ph1*2i)*exp(t*4i) - r^2*z^4*exp(ph1*2i) - 12*r^2*z^3*exp(ph1*1i)*exp(t*2i) + 2*r^2*z^3*exp(ph1*3i)*exp(t*4i) - 4*r^2*z^3*exp(ph1*3i)*exp(t*2i) - 2*r^2*z^3*exp(ph1*1i)*exp(t*4i) + 2*r^2*z^3*exp(ph1*3i) - 2*r^2*z^3*exp(ph1*1i) + 24*r^2*z^2*exp(ph1*2i)*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i)*exp(t*4i) + 2*r^2*z^2*exp(ph1*4i)*exp(t*2i) - 8*base^2*z^2*exp(ph1*2i)*exp(t*2i) - 4*base^2*z^2*exp(ph1*2i)*exp(t*4i) - r^2*z^2*exp(ph1*4i)*exp(t*4i) + 2*r^2*z^2*exp(t*2i) + 4*r^2*z^2*exp(ph1*2i) - 4*base^2*z^2*exp(ph1*2i) - r^2*z^2*exp(t*4i) - r^2*z^2*exp(ph1*4i) - r^2*z^2 - 12*r^2*z*exp(ph1*3i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i)*exp(t*4i) + 2*r^2*z*exp(ph1*1i)*exp(t*4i) - 4*r^2*z*exp(ph1*1i)*exp(t*2i) - 2*r^2*z*exp(ph1*3i) + 2*r^2*z*exp(ph1*1i) + 2*r^2*exp(ph1*2i)*exp(t*2i) - r^2*exp(ph1*2i)*exp(t*4i) - r^2*exp(ph1*2i)==0;
z1 = solve(eq,z,'MaxDegree',4)
%%
matlabFunction(z1,'File','solveQuartic2')