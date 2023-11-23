syms t1t t2t phi3 r x y z
p1t = [-x;y;z];
p2t = [x;y;z];
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
n = Rz*Rt*nt;
center = Rz*Rx*[0 -r 0].'
p1p = p1t-(p1t.'*n)*n;
    %     dist2 = p2t.'*n;
    p2p = p2t-(p2t.'*n)*n;
    %     p12 = p2p-p1p;
    rad1 = p1p-center;
    rad2 = p2p-center;
    eq1 = simplify(rad1.'*rad1-r^2,20)
    eq2 = rad2.'*rad2-r^2
    
    simplify(eq1-eq2)
