function varargout = CylinderCliffPoint(r,hcm,p1,p2,p3,t,phi3,theta12)
x = sym('x',[2,1]);
t1t = x(1);
t2t = x(2);
p2t = p2-p3;
p1t = p1-p3;
% p3t = [0;0;0];
% N = 100;

cliff_height = 0;
cliffTheta1 = 0;
cliffTheta2 = 0;
% n = [sin(t1t)*sin(t2t);-sin(t1t)*cos(t2t);cos(t1t)];
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
n = Rt*nt;
options = optimoptions('fsolve','Display','none');
minh = inf;
    center = Rt*(-r)*[-sin(phi3);cos(phi3);0];
%     center = r*[sin(phi3)*cos(t1t);cos(phi3)*cos(t1t);-sin(t1t)]; %center of the bottom base
%     dist1 = p1t.'*n;
    p1p = p1t-(p1t.'*n)*n;
%     dist2 = p2t.'*n;
    p2p = p2t-(p2t.'*n)*n;
%     p12 = p2p-p1p;
    rad1 = p1p-center;
    rad2 = p2p-center;
    eq1 = rad1.'*rad1-r^2;
    eq2 = rad2.'*rad2-r^2;
        guess = theta12.';

    feq = matlabFunction([eq1;eq2],'vars',{x});
    [sol,~,~,~] = fsolve(feq,guess,options);
    
    comt = p3 + center+n*hcm;
    ht = comt(3);
    cliff_height = subs(ht,[t1t;t2t],sol);
    cliffTheta1 = sol(1);
    cliffTheta2 = sol(2);
    comsol = double(subs(comt,[t1t;t2t],sol));
    if comsol(2)<p2(2)
        fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),p1,p2,p3);
        disp('error in solution')
        close fig;
    end
    if cliff_height<minh
        minCOM = comsol;
        minh = cliff_height;
        minTheta1 =  sol(1);
        minTheta2 =  sol(2);
    end
    
    


%%

switch nargout
    case 0
% save('cliff4.mat','phi3_range','cliff_height','cliffTheta1','cliffTheta2');
figure
hold on
plot(phi3_range,cliff_height)
% plot(cliffTheta1,cliff_height)
xlabel('\phi_3')
ylabel('COM height');
    case 1
        varargout{1} = cliff_height;
    case 2
        varargout{1} = phi3_range;
        varargout{2} = cliff_height;
    case 4
        varargout{1} = phi3_range;
        varargout{2} = cliff_height;
        varargout{3} = cliffTheta1;
        varargout{4} = cliffTheta2;
    case 3
        varargout{1} = minCOM;
        varargout{2} = minTheta1;
        varargout{3} = minTheta2;
end
end

