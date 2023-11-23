syms t p

% contact normals
n1 = -[sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n2 = -[-sin(p),cos(p)*cos(t),cos(p)*sin(t)].';
n3 = [0,-sin(t),cos(t)].';

%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = -N\g

%%
syms ds h3
h = cos(t)*ds;
h1 = h; h2 = h;
h3 = solve(h1*L(1)*n1(2)+h2*L(2)*n2(2)+h3*L(3)*n3(2)==-sin(t)*ds,h3)
simplify(h1*L(1)*n1(1)+h2*L(2)*n2(1))

