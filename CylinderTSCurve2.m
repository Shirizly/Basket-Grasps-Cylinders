function varargout = CylinderTSCurve2(r,d3,hcm,p1,p2,p3,t,phi,N)

% phi3_range = linspace(-pi(),0,N);%linspace(-1.8,-0.8,N);%

options = optimoptions('fsolve','OptimalityTolerance',1E-5,'Display','none','FunctionTolerance',1E-8);

[t1s1,t1s2] = TS_symmetrical2(r,d3,p1,p2,p3,phi); %t1 for phi3=0,pi
phi3_range = linspace(-pi(),0,N);
[circle_height,circleTheta1,circleTheta2,centersol,comsol,nsol] = computeHalfCliffApprox2nd(phi3_range,r,d3,hcm,p1,p2,p3,t,N,t1s2,options);
crit_ind = find(circleTheta2==inf,1,'first');
if ~isempty(crit_ind)
phi3_crit = phi3_range(crit_ind);
circle_height = circle_height(1:crit_ind-1);
circleTheta1 = circleTheta1(1:crit_ind-1);
circleTheta2 = circleTheta2(1:crit_ind-1);
centersol = centersol(1:crit_ind-1,:);
comsol = comsol(1:crit_ind-1,:);
nsol = nsol(1:crit_ind-1,:);
phi3_range = phi3_range(1:crit_ind-1);
phi3_range2 = linspace(0,phi3_crit+pi()/N,N-crit_ind+1);

disp('range 2');
[circle_height2,circleTheta21,circleTheta22,centersol2,comsol2,nsol2] = computeHalfCliffApprox2nd(phi3_range2,r,d3,hcm,p1,p2,p3,t,N-crit_ind+1,-t1s1,options);
crit_ind2 = find(circleTheta22==inf,1,'first');
if ~isempty(crit_ind2)
%     phi3_crit = phi3_range(crit_ind);
% circle_height2 = circle_height2(1:crit_ind2-1);
circleTheta21 = circleTheta21(1:crit_ind2-1);
circleTheta22 = circleTheta22(1:crit_ind2-1);
centersol2 = centersol2(1:crit_ind2-1,:);
comsol2 = comsol2(1:crit_ind2-1,:);
nsol2 = nsol2(1:crit_ind2-1,:);
% phi3_range = phi3_range(1:crit_ind2-1);
end
circle_height = [circle_height circle_height2(end:-1:1)];
circleTheta1 = [circleTheta1 circleTheta21(end:-1:1)];
circleTheta2 = [circleTheta2 circleTheta22(end:-1:1)];
nsol = [nsol;nsol2(end:-1:1,:)];
centersol = [centersol;centersol2(end:-1:1,:)];
comsol = [comsol;comsol2(end:-1:1,:)];
phi3_range = [phi3_range phi3_range2(end:-1:1)];
end

%%
% phi3_range= wrapToPi(phi3_range);
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



function [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2nd(phi3_range,r,d3,p1,p2,p3,t)
cliffHeight = inf*ones(size(phi3_range));
centersolA = inf*ones(size(phi3_range,1),3);
comsolA = centersolA;
nsolA = centersolA;
p2t = p2-p3;
p1t = p1-p3;


[t1c1,t1c2] = TS_symmetrical2(r,d3,p1,p2,p3,phi);
if ~isempty(t1c1)&&~isempty(t1c2)
R1 = [1 0 0; 0 cos(-t1c1) -sin(-t1c1); 0 sin(-t1c1) cos(-t1c1)];
R2 = [1 0 0; 0 cos(-t1c2) -sin(-t1c2); 0 sin(-t1c2) cos(-t1c2)];
t1c1center = R1*[0;-r;0]; %calculated relative to p3
t1c2center = R2*[0;r;0];
circleCenter = (t1c1center+t1c2center)/2;
noon = t1c1center-circleCenter;

% circleAxis = [0;-sin(tmean);cos(tmean)];
% figure
% hold on
perpToBoth = norm(noon)*[-1;0;0];
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
%     Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
    clockArm = noon*cos(phi3) + perpToBoth*sin(phi3);
    center = circleCenter + clockArm; %calculated as if p3==0

    %compute parameters of cones starting from center tangent to spheres
    %surrounding O1,2
    d1 = p1t-center;
    nd1 = norm(d1);
    d1n = d1/nd1;
    ph1 = asin(r/norm(d1));
    
    d2 = p2t-center;
    nd2 = norm(d2);
    d2n = d2/nd2;
    ph2 = asin(r/norm(d2));
    %compute COM position given tangency criterion+base center position
    coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
    sol = coms(2,:); 
%     if ~imag(sol)==0
% %         continue;
%         if imag(coms(1,:))==0
%         sol = coms(1,:);
%         else
%             continue;
%         end
%     end
    pointOnLine = -center.'*sol.'/(sol*sol.')*sol.'+center;
    
%     sol*center
%     norm(pointOnLine)
    center = pointOnLine*r/norm(pointOnLine);
    %some plots to look at certain values relevant to the correction
%     plot(i,r-norm(pointOnLine),'g+')
%     plot(i,sol*center,'r+')
    
%%additional correction  - rerun the computation of the tangent using
    %%the corrected center - turned out to get worse results, sus?
    d1 = p1t-center;
    d1n = d1/norm(d1);
    ph1 = asin(r/norm(d1));
    
    d2 = p2t-center;
    d2n = d2/norm(d2);
    ph2 = asin(r/norm(d2));
    
    coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
    sol = coms(2,:); 
    
    pointOnLine = -center.'*sol.'/(sol*sol.')*sol.'+center;
    center = pointOnLine*r/norm(pointOnLine);



    center = center + p3;
    com = sol.'+center;
    nsol = com-center;
    nsol = nsol/norm(nsol);
    
    nsolA(i,:) = nsol.';
    centersolA(i,:) = center.';
    comsolA(i,:) = com.';
    ht = com(3);
    cliffHeight(i) = ht; 
end
else
    
    
end
if any(imag(cliffHeight)~=0)
    problem_ind = imag(cliffHeight)~=0;
    comsolA(problem_ind,:) = inf;
    centersolA(problem_ind,:) = inf;
    nsolA(problem_ind,:) = inf;
    cliffHeight(problem_ind) = inf;  
end
    
end





function [circleHeight,circleTheta1,circleTheta2,centersolA,comsolA,nsolA] = computeHalfCircle(phi3_range,r,d3,hcm,p1,p2,p3,t,N,t1c,options)
circleHeight = inf*ones(size(phi3_range));
circleTheta1 = zeros(size(phi3_range));
circleTheta2 = inf*ones(size(phi3_range));
x = sym('x',[2,1]);
t1t = x(1);
t2t = x(2);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];

minh = inf;
s3 = norm(p3);
epsi = 1E-2;
h = CylinderComPos(r,hcm,p1,p2,p3,t,0,s3);
h= cos(t);
% supportX = [p1(1), p2(1), p3(1)];
% supportY = [p1(2), p2(2), p3(2)];
centersolA = zeros(N,3);
comsolA = centersolA;
nsolA = centersolA;
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
    Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
    center = Rz*Rx*[0 -d3 0].';
    n = Rz*Rt*nt;
    p1p = p1t-(p1t.'*n)*n;
    p2p = p2t-(p2t.'*n)*n;
    rad1 = p1p-center;
    rad2 = p2p-center;
    eq1 = rad1.'*rad1-r^2;
    eq2 = rad2.'*rad2-r^2;
    
    if i==1
        p12 = (p1t+p2t)/2;
        %extract phi1=phi2 from finger positions and r
        adp = norm(p1t-p2t)/2; %half distance between p1,p2
        bdp = sqrt(r^2-adp^2); %other edge of the right triangle
        phi1 = pi()-atan2(adp,bdp);
        d1 = -r*(cos(phi1));
        dpnorm = norm(p12);
        td30 = -asin(d1/dpnorm)+atan2(-p12(2),p12(3));   
        guess = [t1c*d3/r-td30*(r-d3)/r,0].';
    else
        guess = sol;
    end
    feq = matlabFunction([eq1;eq2],'vars',{x});
    
    [sol,value,flag,~] = fsolve(feq,guess,options);
    
    if flag < 1 || norm(value)>1E-2
%         range = [max(-pi(),sol-max(norm(value)*1000,0.01)),min(pi(),sol+max(norm(value)*1000,0.01))];
%         [fig,minval,minpoint] = evalSurface(feq,range);
%         normval = norm(value)
%         minval
%         [sol,value,flag,~] = fsolve(feq,minpoint,options);
%         disp('found critical point, splitting range');
%         close(fig)
%         if flag < 1 || norm(value)>1E-5
            return
%         end
    end
    
    

    Rxsol = [1 0 0; 0 cos(sol(1)) -sin(sol(1)); 0 sin(sol(1)) cos(sol(1))];
    Rysol = [cos(sol(2)) 0 sin(sol(2));0 1 0; -sin(sol(2)) 0 cos(sol(2))];
    Rtsol = Rxsol*Rysol;
    nsol = Rz*Rtsol*nt;
    nsolA(i,:) = nsol.';
    centersol = Rz*Rxsol*[0;-d3;0];
    centersolA(i,:) = centersol.';
    comsol = p3 + centersol+nsol*hcm;
    comsolA(i,:) = comsol.';
    ht = comsol(3);
    circleHeight(i) = ht;
    circleTheta1(i) = sol(1);
    circleTheta2(i) = sol(2);
    %         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);

    dist1 = (p1t-centersol).'*nsol;
    dist2 = (p2t-centersol).'*nsol;
    if dist1*dist2<0
        circleHeight(i)=inf;
    end
    %% check if normals, projected on x-y plane, define a convex cone that includes the origin 
%     [Normals,L,W] = TSNormals(r,hcm,comsol,phi3,p1,p2,p3,sol(1),sol(2));
% 
%     A = Normals(1:2,:);
%     NormalZ = abs(Normals(3,:));
%     rA = rank(A.'*A);
%     if rA<2
%         disp('suspicious rank of force magnitudes')
%     end
%     nullA = null(A);
%     if any(nullA<0)&&any(nullA>0)
% %          circleHeight(i) = inf;
%     end  

%   fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
    if circleHeight(i)<h
         disp('unstable')
    end
end
end


