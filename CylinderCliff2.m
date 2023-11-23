function varargout = CylinderCliff2(r,hcm,phi,p1,p2,p3,t,N)

% p3t = [0;0;0];
% N = 100;
epsi = 1E-2;
options = optimoptions('fsolve','OptimalityTolerance',1E-5,'Display','none','FunctionTolerance',1E-8);
[t1c1,t1c2] = TS_symmetrical(r,p1,p2,p3,phi);

% tic
% optionsf = optimoptions('fmincon','OptimalityTolerance',1E-5,'Display','none','FunctionTolerance',1E-8);
% [cliffPoint,flag] = computeCliffCriticalPointsFmin(r,hcm,p1,p2,p3,t,[t1c1,t1c2],optionsf)
% toc
t1range = [t1c1,t1c2];
phi3_range = linspace(-pi(),0,N);
phi3_rangeA = phi3_range;
tic
if norm(p1-p3)>2*r
[cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2nd(phi3_rangeA,r,p1,p2,p3,t,phi);
else
    [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndSpecial(phi3_rangeA,r,p1,p2,p3,t);
end
toc

% phi3_range = linspace(0,-pi(),N);
tic
[cliff_height,cliffTheta1,cliffTheta2,centersol,comsol,nsol] = computeHalfCliff(phi3_range,r,hcm,p1,p2,p3,t,N,t1c2,options);
toc
crit_ind = find(cliffTheta2==inf,1,'first');
if ~isempty(crit_ind)
phi3_crit = phi3_range(crit_ind);
cliff_height = cliff_height(1:crit_ind-1);
cliffTheta1 = cliffTheta1(1:crit_ind-1);
cliffTheta2 = cliffTheta2(1:crit_ind-1);
centersol = centersol(1:crit_ind-1,:);
comsol = comsol(1:crit_ind-1,:);
nsol = nsol(1:crit_ind-1,:);
phi3_range = phi3_range(1:crit_ind-1);
phi3_range2 = linspace(0,phi3_crit+pi()/N,N-crit_ind+1);

disp('range 2');
[cliff_height2,cliffTheta21,cliffTheta22,centersol2,comsol2,nsol2] = computeHalfCliff(phi3_range2,r,hcm,p1,p2,p3,t,N-crit_ind+1,-t1c1,options);
crit_ind2 = find(cliffTheta22==inf,1,'first');
% if ~isempty(crit_ind2)
% %     phi3_crit = phi3_range(crit_ind);
% cliff_height2 = cliff_height2(1:crit_ind2-1);
% cliffTheta21 = cliffTheta21(1:crit_ind2-1);
% cliffTheta22 = cliffTheta22(1:crit_ind2-1);
% centersol2 = centersol2(1:crit_ind2-1,:);
% comsol2 = comsol2(1:crit_ind2-1,:);
% nsol2 = nsol2(1:crit_ind2-1,:);
% phi3_range2 = phi3_range2(1:crit_ind2-1);
% end


cliff_height = [cliff_height cliff_height2(end:-1:1)];
cliffTheta1 = [cliffTheta1 cliffTheta21(end:-1:1)];
cliffTheta2 = [cliffTheta2 cliffTheta22(end:-1:1)];
nsol = [nsol;nsol2(end:-1:1,:)];
centersol = [centersol;centersol2(end:-1:1,:)];
comsol = [comsol;comsol2(end:-1:1,:)];
phi3_range = [phi3_range phi3_range2(end:-1:1)];

end
toc
[~,minind] = min(cliff_height);
minCOM = comsol(minind,:).';
minTheta1 = cliffTheta1(minind);
minTheta2 = cliffTheta2(minind);
minPhi3 = phi3_range(minind);
drawq = 0;
if minind>1&&minind<N
    [minTheta1;minTheta2]
    guess = [minTheta1,minTheta2,minPhi3,1,1];
    cliffPoints = computeCliffCriticalPoints(r,hcm,p1,p2,p3,t,guess,options)
    disp('non-symmetrical');
    drawq = 1;
end
difference = zeros(N,1);
comsol(:,1) = -comsol(:,1);
for i=1:N
    p1A = comsolA(i,:);
    dist = comsol-p1A;
    dist_sq = diag(dist*dist.');
    ind = find(dist_sq==min(dist_sq));
    p1c = comsol(ind,:);
    if ind==1 || ind == size(comsol,1)
        difference(i) = norm(p1A-p1c);
    else
        p1Arel = p1A-p1c;
        p21c = comsol(ind-1,:);
        p22c = comsol(ind+1,:);
        v1 = p21c-p1c;
        v1 = v1/norm(v1);
        v2 = p22c-p1;
        v2 = v2/norm(v2);
        
        dist1 = norm((eye(3)-v1.'*v1)*p1Arel.');
        dist2 = norm((eye(3)-v2.'*v2)*p1Arel.');
        difference(i) = min(dist1,dist2);
    end
        
end

if 1

figure
hold on
plot3(centersol(:,1),centersol(:,2),centersol(:,3),'-*')
% plot3(centersol2(:,1),centersol2(:,2),centersol2(:,3),'-*')
xlabel('x')
ylabel('y')
zlabel('z')

figure
hold on
plot3(comsol(:,1),comsol(:,2),comsol(:,3),'-*')
% plot3(comsol2(:,1),comsol2(:,2),comsol2(:,3),'-*')
xlabel('x')
ylabel('y')
zlabel('z')
end

difference_mean = mean(difference)
% if drawq || difference_mean>1E-1
f = figure;
    hold on
    plot3(p1(1),p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
    plot3(p2(1),p2(2),p2(3),'ob','markerSize',4,'lineWidth',4);
    plot3(p3(1),p3(2),p3(3),'og','markerSize',4,'lineWidth',4);
    [x,y,z] = sphere;
    mesh(r*x+p1(1),r*y+p1(2),r*z+p1(3),'faceAlpha',0.2);
    mesh(r*x+p2(1),r*y+p2(2),r*z+p2(3),'faceAlpha',0.2);
%     mesh(r*x+p3(1),r*y+p3(2),r*z+p3(3),'faceAlpha',0.2);
c = centersol(20,:).';
[X,Y,Z]=cylinder([0 .5],50 );
axis([0 1,-1 1,-.5 .5])
M=makehgtform('translate',[0,0,0],'xrotate',pi/4,'yrotate',pi/2);
h=surf(X+c(1),Y+c(2),Z+c(3),'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.4);

grid on



f = figure;
    hold on
    plot3(p1(1),p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
    plot3(p2(1),p2(2),p2(3),'ob','markerSize',4,'lineWidth',4);
    plot3(p3(1),p3(2),p3(3),'og','markerSize',4,'lineWidth',4);
    [x,y,z] = sphere;
    mesh(r*x+p1(1),r*y+p1(2),r*z+p1(3),'faceAlpha',0.2);
    mesh(r*x+p2(1),r*y+p2(2),r*z+p2(3),'faceAlpha',0.2);
    mesh(r*x+p3(1),r*y+p3(2),r*z+p3(3),'faceAlpha',0.2);
    
    for i=1:size(nsol,1)
        ax = nsol(i,:).';
        center = p3 + centersol(i,:).';
        axisends = [center center+2*ax];
        plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k','lineWidth',1);
        plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k','lineWidth',4,'markerSize',2);
        axisends = center+ax;
        plot3(axisends(1,:),axisends(2,:),axisends(3,:),'ro','lineWidth',5,'markerSize',2);
        ax(1) = -ax(1);
        center(1) = -center(1);
        axisends = [center center+2*ax];
        plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k','lineWidth',1);
        plot3(axisends(1,:),axisends(2,:),axisends(3,:),'k','lineWidth',4,'markerSize',2);
        axisends = center+ax;
        plot3(axisends(1,:),axisends(2,:),axisends(3,:),'ro','lineWidth',5,'markerSize',2);
        
        centerA = centersolA(i,:).';
        comA = comsolA(i,:).';
        axisA = [centerA comA];
        
%         plot3(axisA(1,:),axisA(2,:),axisA(3,:),'ro-','lineWidth',4,'markerSize',2);
    end
    
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')


% if difference_mean>1E-1
    f2 = figure;
    hold on
    plot(phi3_range,cliff_height,'b')
    plot(phi3_rangeA,cliffHeight,'k')
    xlabel('\phi_3')
    ylabel('COM height');
    
    % angle1 = zeros(length(cliff_height),1);
    % for i=1:length(cliff_height)
    %     angle1(i) = atan2(comsol(i,2),comsol(i,1));
    % end
    % angle2 = zeros(length(cliffHeight),1);
    % for i=1:length(cliffHeight)
    %     angle2(i) = atan2(comsolA(i,2),-comsolA(i,1));
    % end
    % f3 = figure;
    % hold on
    % plot(angle1,cliff_height,'b')
    % plot(angle2,cliffHeight,'k')
    % xlabel('angle')
    % ylabel('COM height');
    
    f4 = figure;
    hold on
    plot3(comsol(:,1),comsol(:,2),comsol(:,3),'b')
    plot3(comsolA(:,1),comsolA(:,2),comsolA(:,3),'k+')
    legend('Accurate','Approximated')
    xlabel('x')
    ylabel('y')
    zlabel('z')
% end
%%
% guess = [minTheta1,minTheta2,minPhi3,1,1];
% [cliffPoints,flag] = computeCliffCriticalPoints(r,hcm,p1,p2,p3,t,guess,options);

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
        fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,minCOM,minTheta1,minTheta2,minPhi3,p1,p2,p3);
    case 1
        varargout{1} = cliff_height;
    case 2
        varargout{1} = phi3_range;
        varargout{2} = cliff_height;
    case 4
        varargout{1} = phi3_range;
        varargout{2} = cliff_height;
        varargout{3} = cliffHeight;%cliffTheta1;
        varargout{4} = cliffTheta2;
    case 3
        varargout{1} = minCOM;
        varargout{2} = minTheta1;
        varargout{3} = minTheta2;
end
end


function [cliffPoint,flag] = computeCliffCriticalPoints(r,hcm,p1,p2,p3,t,guess,options)
x = sym('x',[3,1]);
syms lam1 lam2
t1t = x(1);
t2t = x(2);
phi3 = x(3);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rz*Rx*[0 -r 0].';
n = Rz*Rt*nt;
p1p = p1t-(p1t.'*n)*n;
p2p = p2t-(p2t.'*n)*n;
rad1 = p1p-center;
rad2 = p2p-center;
eq1 = rad1.'*rad1-r^2;
eq2 = rad2.'*rad2-r^2;
com = p3 + center+n*hcm;
comgrad = gradient(com(3),x());
eq1grad = gradient(eq1,x);
eq2grad = gradient(eq2,x);
xs = [x;lam1;lam2];
feq = matlabFunction([eq1;eq2;lam1*eq1grad+lam2*eq2grad-comgrad],'vars',{xs});

[sol,~,flag,~] = fsolve(feq,guess.',options);

Rxsol = [1 0 0; 0 cos(sol(1)) -sin(sol(1)); 0 sin(sol(1)) cos(sol(1))];
Rysol = [cos(sol(2)) 0 sin(sol(2));0 1 0; -sin(sol(2)) 0 cos(sol(2))];
Rtsol = Rxsol*Rysol;
Rz = [cos(sol(3)) -sin(sol(3)) 0;sin(sol(3)) cos(sol(3)) 0;0 0 1];
nsol = Rz*Rtsol*nt;
centersol = Rz*Rxsol*[0;-r;0];
comsol = p3 + centersol+nsol*hcm;
cliffPoint = [sol;comsol];
end

function [cliffPoint,flag] = computeCliffCriticalPointsFmin(r,hcm,p1,p2,p3,t,tc,options)
x = sym('x',[3,1]);

t1t = x(1);
t2t = x(2);
phi3 = x(3);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rz*Rx*[0 -r 0].';
n = Rz*Rt*nt;
% the non-linear conditions for the cliff curve.
p1p = p1t-(p1t.'*n)*n;
p2p = p2t-(p2t.'*n)*n;
rad1 = p1p-center;
rad2 = p2p-center;
eq1 = rad1.'*rad1-r^2; %equation for contact of finger 1 with cyl. surf.
eq2 = rad2.'*rad2-r^2;
dist1 = (p1t-center).'*n; %distance of finger 1 from the base of the cyl.
dist2 = (p2t-center).'*n;
p12 = (p1t+p2t)/2;
con3 = cross(center-p12,n);
nonlcon = matlabFunction([-dist1;-dist2;-con3(1)],[eq1;eq2],'vars',{x},'Outputs',{'c','ceq'});

com = p3 + center+n*hcm;
fun = matlabFunction(com(3),'vars',{x});
guess1 = [-tc(1),0,0];
[solt(:,1),valuet(1),flag,output] = fmincon(fun,guess1.',[],[],[],[],[],[],nonlcon,options);
guess2 = [tc(2),0,0];
[solt(:,2),valuet(2),flag2,output2] = fmincon(fun,guess2.',[],[],[],[],[],[],nonlcon,options);
guess3 = [t,0,0];
[solt(:,3),valuet(3),flag2,output2] = fmincon(fun,guess3.',[],[],[],[],[],[],nonlcon,options);
[value,min_ind] =min(valuet);
sol = solt(:,min_ind);


Rxsol = [1 0 0; 0 cos(sol(1)) -sin(sol(1)); 0 sin(sol(1)) cos(sol(1))];
Rysol = [cos(sol(2)) 0 sin(sol(2));0 1 0; -sin(sol(2)) 0 cos(sol(2))];
Rtsol = Rxsol*Rysol;
Rz = [cos(sol(3)) -sin(sol(3)) 0;sin(sol(3)) cos(sol(3)) 0;0 0 1];
nsol = Rz*Rtsol*nt;
centersol = Rz*Rxsol*[0;-r;0];
comsol = p3 + centersol+nsol*hcm;
dist1 = (p1t-centersol).'*nsol;
    dist2 = (p2t-centersol).'*nsol;
cliffPoint = [sol;comsol];
end

function [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox(phi3_range,r,p1,p2,p3,t)
cliffHeight = zeros(size(phi3_range));
centersolA = zeros(size(phi3_range,1),3);
comsolA = centersolA;
nsolA = centersolA;
p2t = p2-p3;
p1t = p1-p3;


[t1c1,t1c2] = TS_symmetrical(r,p1,p2,p3);
tmean = (t1c1+t1c2)/2;
R1 = [1 0 0; 0 cos(-t1c1) -sin(-t1c1); 0 sin(-t1c1) cos(-t1c1)];
R2 = [1 0 0; 0 cos(-t1c2) -sin(-t1c2); 0 sin(-t1c2) cos(-t1c2)];
t1c1center = R1*[0;-r;0]; %calculated relative to p3
t1c2center = R2*[0;r;0];
circleCenter = (t1c1center+t1c2center)/2;
noon = t1c1center-circleCenter;

% circleAxis = [0;-sin(tmean);cos(tmean)];

perpToBoth = norm(noon)*[-1;0;0];
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
%     Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
    clockArm = noon*cos(phi3) + perpToBoth*sin(phi3);
    center = circleCenter + clockArm;
    
    d1 = p1t-center;
    d1n = d1/norm(d1);
    ph1 = asin(r/norm(d1));
    
    d2 = p2t-center;
    d2n = d2/norm(d2);
    ph2 = asin(r/norm(d2));
    
    coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
    sol = coms(2,:); 
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
end

function [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2nd(phi3_range,r,p1,p2,p3,t,phi)
cliffHeight = inf*ones(size(phi3_range));
centersolA = inf*ones(size(phi3_range,1),3);
comsolA = centersolA;
nsolA = centersolA;
p2t = p2-p3;
p1t = p1-p3;


[t1c1,t1c2] = TS_symmetrical(r,p1,p2,p3,phi);
tmean = (t1c1+t1c2)/2;
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
tic
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
%     d1 = p1t-center;
%     d1n = d1/norm(d1);
%     ph1 = asin(r/norm(d1));
%     
%     d2 = p2t-center;
%     d2n = d2/norm(d2);
%     ph2 = asin(r/norm(d2));
%     
%     coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
%     sol = coms(2,:); 
%     
%     pointOnLine = -center.'*sol.'/(sol*sol.')*sol.'+center;
%     center = pointOnLine*r/norm(pointOnLine);



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
if any(imag(cliffHeight)~=0)
    problem_ind = imag(cliffHeight)~=0;
    comsolA(problem_ind,:) = inf;
    centersolA(problem_ind,:) = inf;
    nsolA(problem_ind,:) = inf;
    cliffHeight(problem_ind) = inf;  
end
toc    
end

function [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndSpecial(phi3_range,r,p1,p2,p3,t)
cliffHeight = inf*ones(size(phi3_range));
centersolA = inf*ones(size(phi3_range,1),3);
comsolA = centersolA;
nsolA = centersolA;
p2t = p2-p3;
p1t = p1-p3;


[t1c1,t1c2] = TS_symmetrical(r,p1,p2,p3);
tmean = (t1c1+t1c2)/2;
R1 = [1 0 0; 0 cos(-t1c1) -sin(-t1c1); 0 sin(-t1c1) cos(-t1c1)];
R2 = [1 0 0; 0 cos(-t1c2) -sin(-t1c2); 0 sin(-t1c2) cos(-t1c2)];
t1c1center = R1*[0;-r;0]; %calculated relative to p3
t1c2center = R2*[0;r;0];
circleCenter = (t1c1center+t1c2center)/2;
noon = t1c1center-circleCenter;

% circleAxis = [0;-sin(tmean);cos(tmean)];

perpToBoth = norm(noon)*[-1;0;0];
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
%     Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
    clockArm = noon*cos(phi3) + perpToBoth*sin(phi3);
    center = circleCenter + clockArm; %calculated as if p3==0
    % check if the base center is within the sphere of one of the other
    % fingers: if yes, skip the calculation
    d1 = p1t-center;
    nd1 = norm(d1);
    d2 = p2t-center;
    nd2 = norm(d2);
    if nd1<r || nd2<r
        continue;
    end
    
    %compute parameters of cones starting from center tangent to spheres
    %surrounding O1,2
    d1n = d1/nd1;
    ph1 = asin(r/norm(d1));
    
    d2n = d2/nd2;
    ph2 = asin(r/norm(d2));
    %compute COM position given tangency criterion+base center position
    coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
    sol = coms(2,:); 
    if ~imag(sol)==0
%         continue;
        if imag(coms(1,:))==0
        sol = coms(1,:);
        else
            continue;
        end
    end
    pointOnLine = -center.'*sol.'/(sol*sol.')*sol.'+center;
%     sol*center
%     norm(pointOnLine)
    center = pointOnLine*r/norm(pointOnLine);
    
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
%     norm(pointOnLine)
%     sol*center
    if norm(p1t-center)<r || norm(p2t-center)<r
        continue;
    end

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
if any(imag(cliffHeight)~=0)
    problem_ind = imag(cliffHeight)~=0;
    comsolA(problem_ind,:) = inf;
    centersolA(problem_ind,:) = inf;
    nsolA(problem_ind,:) = inf;
    cliffHeight(problem_ind) = inf;  
end
    
end


function [cliffHeight,cliffTheta1,cliffTheta2,centersolA,comsolA,nsolA] = computeHalfCliff(phi3_range,r,hcm,p1,p2,p3,t,N,t1c,options)

cliffHeight = zeros(size(phi3_range));
cliffTheta1 = zeros(size(phi3_range));
cliffTheta2 = inf*ones(size(phi3_range));
x = sym('x',[2,1]);
t1t = x(1);
t2t = x(2);
p2t = p2-p3;
p1t = p1-p3;
% n = [sin(t1t)*sin(t2t);-sin(t1t)*cos(t2t);cos(t1t)];
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];

minh = inf;
s3 = norm(p3);
epsi = 1E-2;
P = [p1,p2,p3];
h = CylinderComPos(r,hcm,P,t,0,s3)
% supportX = [p1(1), p2(1), p3(1)];
% supportY = [p1(2), p2(2), p3(2)];
centersolA = zeros(N,3);
comsolA = centersolA;
nsolA = centersolA;
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
    Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
    center = Rz*Rx*[0 -r 0].';
%     center = (-r)*[-sin(phi3)*cos(t1t),cos(t1t)*cos(phi3),-sin(t1t)].';
    n = Rz*Rt*nt;
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
        guess = [t1c,0].';
    else
        guess = sol;
    end
    feq = matlabFunction([eq1;eq2],'vars',{x});
    
    [sol,value,flag,~] = fsolve(feq,guess,options);
    
    while flag < 1 || norm(value)>1E-4
%         range = [max(-pi(),sol-max(norm(value)*1000,0.01)),min(pi(),sol+max(norm(value)*1000,0.01))];
%         [fig,minval,minpoint] = evalSurface(feq,range);
%         normval = norm(value)
%         minval
%         [sol,value,flag,~] = fsolve(feq,minpoint,options);
%         disp('found critical point, splitting range');
%         close(fig)
%         if flag<1
            return
%         end
    end
    
    
    % this check is wrong due to different definitions of theta1, probably
    % can be fixed
%     if (sol(1)>(max(t1range)+epsi))||sol(1)<(min(t1range)-epsi)
%         disp('probable error')
%     end
    Rxsol = [1 0 0; 0 cos(sol(1)) -sin(sol(1)); 0 sin(sol(1)) cos(sol(1))];
    Rysol = [cos(sol(2)) 0 sin(sol(2));0 1 0; -sin(sol(2)) 0 cos(sol(2))];
    Rtsol = Rxsol*Rysol;
    nsol = Rz*Rtsol*nt;
    nsolA(i,:) = nsol.';
    centersol = Rz*Rxsol*[0;-r;0];
    centersolA(i,:) = centersol.';
    comsol = p3 + centersol+nsol*hcm;
    comsolA(i,:) = comsol.';
    ht = comsol(3);
    cliffHeight(i) = ht;
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
%         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
%         disp('error in solution')
%         close(fig)
%     end
    %% check if one of the fingers 1,2 reached the edge of the bottom base. If yes, everything beyond has height inf (non-feasible).

    dist1 = (p1t-centersol).'*nsol;
    dist2 = (p2t-centersol).'*nsol;
    if dist1*dist2<0
        cliffHeight(i)=inf;
    end
    %% check if normals, projected on x-y plane, define a convex cone that includes the origin 
%         [Normals,L,W] = TSNormals(r,hcm,comsol,phi3,p1,p2,p3,sol(1),sol(2));

%     A = Normals(1:2,:);
%     rA = rank(A.'*A);
%     if rA<2
%         disp('suspicious rank of force magnitudes')
%     end
%     nullA = null(A);
%     if any(nullA<0)&&any(nullA>0)
%          cliffHeight(i) = inf;
%     end
    
    
%     if 
%         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
%         disp('error in solution')
%     end
    if cliffHeight(i)<minh
        minCOM = comsol;
        minh = cliffHeight(i);
        minTheta1 =  sol(1);
        minTheta2 =  sol(2);
        minPhi3 = phi3;
    end
    
    if minh<h
        disp('unstable')
    end
end

end



function [fig,minval,minpoint] = evalSurface(f,range)
% x = -pi():0.03:pi();
% y = -pi():0.03:pi();
x = linspace(range(1,1),range(1,2),201);
y = linspace(range(2,1),range(2,2),201);
fxy = zeros(length(x),length(y));
for i=1:length(x)
    xt = x(i);
    for j=1:length(y)
        yt = y(j);
        fxyt = feval(f,[xt;yt]);
        fxy(i,j) = norm(fxyt);
    end
end
fig = figure
surf(x,y,fxy.')
minval = min(min(fxy));
minpoint_ind = find(fxy==minval);
y_ind = mod(minpoint_ind,length(x));
x_ind = (minpoint_ind-y_ind)/length(x)+1;
minpoint = [x(x_ind);y(y_ind)];
end
        