
phi3 = phi3_range(minind);
    center = Rt*(-r)*[-sin(phi3);cos(phi3);0]; % center of the bottom base relative to p3
    %     dist1 = p1t.'*n;
    p1p = p1t-(p1t.'*n)*n;
    %     dist2 = p2t.'*n;
    p2p = p2t-(p2t.'*n)*n;
    %     p12 = p2p-p1p;
    rad1 = p1p-center;
    rad2 = p2p-center;
    eq1 = rad1.'*rad1-r^2;
    eq2 = rad2.'*rad2-r^2;
    if i==1 %|| dist1*dist2<0
        guess = [t,1E-5].';
    else
        guess = sol;
    end
    feq = matlabFunction([eq1;eq2],'vars',{x});
    [sol,~,~,~] = fsolve(feq,guess,options);
    Rxsol = [1 0 0; 0 cos(sol(1)) -sin(sol(1)); 0 sin(sol(1)) cos(sol(1))];
    Rysol = [cos(sol(2)) 0 sin(sol(2));0 1 0; -sin(sol(2)) 0 cos(sol(2))];
    Rtsol = Rxsol*Rysol;
    nsol = Rtsol*nt;
    centersol = Rtsol*(-r)*[-sin(phi3);cos(phi3);0];
    comsol = p3 + centersol+nsol*hcm;
    ht = comsol(3);

        fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
