function varargout = CylinderTSCurve3(r,d3,hcm,P,t,phi,N,prev_sol,feq2,t2_sol)

% phi3_range = linspace(-pi(),0,N);%linspace(-1.8,-0.8,N);%


options = optimoptions('fsolve','OptimalityTolerance',1E-5,'Display','none','FunctionTolerance',1E-8);

% [t1s1,t1s2] = TS_symmetrical2(r,d3,P,phi); %t1 for phi3=0,pi
phi3_range = linspace(-pi(),0,N);
[circle_height,circleTheta1,circleTheta2,centersol,comsol,nsol] = computeHalfCircle(phi3_range,r,d3,hcm,P,t,N,options,prev_sol,feq2,t2_sol);


%%

switch nargout
    case 0
        % save('circle4.mat','phi3_range','circle_height','circleTheta1','circleTheta2');
        figure
        hold on
        plot(phi3_range,circle_height)
        % plot(circleTheta1,circle_height)
        xlabel('\phi_3')
        ylabel('COM height');
        fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,minCOM,minTheta1,minTheta2,p1,p2,p3);
    case 1
        varargout{1} = circle_height;
    case 2
        varargout{1} = phi3_range;
        varargout{2} = circle_height;
    case 4
        varargout{1} = phi3_range;
        varargout{2} = circle_height;
        varargout{3} = circleTheta1;
        varargout{4} = circleTheta2;
    case 3
        varargout{1} = minCOM;
        varargout{2} = minTheta1;
        varargout{3} = minTheta2;
end
end

function [circleHeight,circleTheta1,circleTheta2,centersolA,comsolA,nsolA] = computeHalfCircle(phi3_range,r,d3,hcm,P,t,N,options,prev_sol,feq2,t2_sol)
circleHeight = inf*ones(size(phi3_range));
circleTheta1 = zeros(size(phi3_range));
circleTheta2 = inf*ones(size(phi3_range));
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
minh = inf;
s3 = norm(p3);
epsi = 1E-2;
h = CylinderComPos(r,hcm,P,t,0,s3);
h= cos(t);
% supportX = [p1(1), p2(1), p3(1)];
% supportY = [p1(2), p2(2), p3(2)];
centersolA = zeros(N,3);
comsolA = centersolA;
nsolA = centersolA;
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
    feq1 = @(x) feq2(x,phi3);
    guess = prev_sol(:,i); %take as initial guess the solution from the same angle, one step close to the origin
    if any(guess==inf)
        guess = [t;t];
    end

    try
    [sol,value,flag,~] = fsolve(feq1,guess,options);
    
    catch
        disp('error in fsolve')
    end
    if flag < 1 || norm(value)>1E-2
        sol = [];
        continue
    end
    
 for j=1:2   

    Rxsol = [1 0 0; 0 cos(sol(1)) -sin(sol(1)); 0 sin(sol(1)) cos(sol(1))];
    Rysol = [cos(sol(2)) 0 sin(sol(2));0 1 0; -sin(sol(2)) 0 cos(sol(2))];
    Rtsol = Rxsol*Rysol;
    nsol = Rtsol*[0;0;1];
    nsolA(i,:) = nsol.';
    centersol = p3 + Rxsol*Rysol*[d3*sin(phi3);-d3*cos(phi3);0];
    centersolA(i,:) = centersol.';
    comsol = centersol+nsol*hcm;
    comsolA(i,:) = comsol.';
    ht = comsol(3);
    circleHeight(i) = ht;
    circleTheta1(i) = sol(1);
    circleTheta2(i) = sol(2);
    %         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
    %% check if center of mass outside the triangle defined by the fingers
    if 0
    xq = comsol(1);
    yq = comsol(2);
    xv = [p1(1) p2(1) p3(1)];
    yv = [p1(2) p2(2) p3(2)];
    if ~inpolygon(xq,yq,xv,yv)
        circleHeight(i) = inf;
        [N,L,W] = TSNormals(r,hcm,comsol,phi3,p1,p2,p3,sol(1),sol(2));
        if all(L>0)
            disp('feasible eq in expected infeasible eq region');
        end
    end
    end
    
    %% check if the side fingers crossed the bases edge
    dist1 = (p1-centersol).'*nsol;
    dist2 = (p2-centersol).'*nsol;
    if any([dist1,dist2]<0)||any([dist1,dist2]>2*hcm)
        circleHeight(i)=inf;
    end
%     if any([dist1,dist2]>2*hcm)
%         disp('fingers too high')
%     end
    %% check if axis in solution "crossed over" point between fingers
    % if it did, that implies that the relevant solution doesn't exist.
    p12 = (p1+p2)/2;
    comsolp = comsol;
    comsolp(1) = 0;
    centersolp = centersol;
    centersolp(1) = 0;
    nsolp = nsol;
    nsolp(1) = 0;
    s = (comsolp-centersolp).'*(centersolp-p12)/norm(comsolp-centersolp)^2;
    arm = (centersolp-s*(comsolp-centersolp))-p12;
    torque = cross(arm,nsolp);
    if torque(1)<0
        circleHeight(i) = inf;
%         [dist1,dist2]
    end
    x = d3*sin(phi3);
    y = d3*cos(phi3);

%   fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
    if circleHeight(i)<h
         disp('unstable')
    end
    if circleHeight(i) ~= inf
        break
    end
end
end
end


