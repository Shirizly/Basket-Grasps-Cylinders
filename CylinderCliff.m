function varargout = CylinderCliff(r,hcm,p1,p2,p3,t,N)
x = sym('x',[2,1]);
t1t = x(1);
t2t = x(2);
p2t = p2-p3;
p1t = p1-p3;
% p3t = [0;0;0];
% N = 100;
phi3_range = linspace(-pi(),0,N);
cliff_height = zeros(size(phi3_range));
cliffTheta1 = zeros(size(phi3_range));
cliffTheta2 = zeros(size(phi3_range));
% n = [sin(t1t)*sin(t2t);-sin(t1t)*cos(t2t);cos(t1t)];
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
n = Rt*nt;
options = optimoptions('fsolve','Display','none');
minh = inf;
s3 = norm(p3);
h = CylinderComPos(r,hcm,p1,p2,p3,t,0,s3);
supportX = [p1(1), p2(1), p3(1)];
supportY = [p1(2), p2(2), p3(2)];
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
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
    cliff_height(i) = ht;
    cliffTheta1(i) = sol(1);
    cliffTheta2(i) = sol(2);
    %     if ~inpolygon(comsol(1),comsol(2),supportX,supportY)
    %         cliff_height(i) = inf;
    %     end
%     if comsol(2)*p2(2)<0
%         cliff_height(i) = inf;
%     end
%     if cliffTheta1(i)>pi()/2 || cliffTheta1(i)<0 %then the cliff is skipping through a single or double-support escape stance
%         cliff_height(i) = inf;
%     end
%     if comsol(2)<p2(2) % the center of mass should never reach "behind" the fingers, unless the solver found the infeasible solution (or rare case)
%         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),p1,p2,p3);
%         disp('error in solution')
%         close(fig)
%     end
    %% check if one of the fingers 1,2 reached the edge of the bottom base. If yes, everything beyond has height inf (non-feasible).

    dist1 = (p1t-centersol).'*nsol;
    dist2 = (p2t-centersol).'*nsol;
    if dist1*dist2<0
        cliff_height(i)=inf;
    end
%     if abs(sol(2))>pi()/4
%         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),p1,p2,p3);
%         disp('error in solution')
%     end
    if cliff_height(i)<minh 
        minCOM = comsol;
        minh = cliff_height(i);
        minTheta1 =  sol(1);
        minTheta2 =  sol(2);
    end
    
    if minh<h
        disp('unstable')
    end
end
[~,minind] = min(cliff_height);
if minind>1&&minind<N
    disp('non-symmetrical');
end
%%
% phi3_range= wrapToPi(phi3_range);
switch nargout
    case 0
        % save('cliff4.mat','phi3_range','cliff_height','cliffTheta1','cliffTheta2');
        figure
        hold on
        plot(phi3_range,cliff_height)
        % plot(cliffTheta1,cliff_height)
        xlabel('\phi_3')
        ylabel('COM height');
        fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,minCOM,minTheta1,minTheta2,p1,p2,p3);
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

