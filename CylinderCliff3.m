function varargout = CylinderCliff3(r,hcm,p1,p2,p3,t,phi,N)
% This function returns the cliff curve, a curve of stances that forms the
% 3-contact strata of the rim curve of the basket grasp. Any point along
% the cliff curve allows escape from the basin. 
% The cliff curve may consist of several different parts:
% a) the base finger reaching the base edge
% b) one of the surface fingers reaching the bottom base
% c) one of the surface fingers reaching the top base
draw = 0;
% N = 100;
epsi = 1E-4;
% options = optimoptions('fsolve','OptimalityTolerance',1E-5,'Display','none','FunctionTolerance',1E-8);

% first, compute the symmetric cliff stances, assuming they exist
[t1s1,t1s2] = TS_symmetrical(r,p1,p2,p3,phi);
%% compute the cliff curve
phi3_range = linspace(-pi(),pi(),N);
phi3_rangeA = phi3_range;
iterate = 1;
[cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC3(phi3_rangeA,r,p1,p2,p3,t1s1,t1s2,iterate);
for fake_counter =  []
% hopefully complete the approximated solutions for other cases, for now
% only the easy case
%% Identify the type of cliff curve to be computed - currently only the complete circular kind is implemented
cond_list = zeros(4,1);

if ~isempty(t1s1) %check that the first stance exists
   cond_list(1) = 1;
end
dpnorm = norm(0.5*(p1+p2)-p3);
cond_list(2) = 2*((r*(1+cos(phi)))^2+(2*hcm)^2-dpnorm^2>=0); % check that second stance exists - that the cylinder is tall enough




switch sum(cond_list(1:2))
    case 3 %all conditions fulfilled - cliff curve from t1s1 to t1s2
        [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC3(phi3_rangeA,r,p1,p2,p3,t1s1,t1s2);
    case 1 %first condition fulfilled:
        %find a midpoint along the curve where side finger reache base edge
        midpoints = TS_sidefingeredge(r,p1,p2,p3,phi);
        if isempty(midpoints)
            cliffHeight = inf;
            centersolA = 0;
            comsolA = 0;
            nsolA = 0;
        else
        % naive cliff curve from t1s1 to midpoint:
        [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC2(phi3_rangeA,r,p1,p2,p3,t1s1,midpoint);
        % followed by special curve 1 from midpoint to symm. endpoint:
        end
    case 2 %second condition fulfilled
        midpoints = TS_sidefingeredge(r,p1,p2,p3,phi);
        if isempty(midpoints)
            cliffHeight = inf;
            centersolA = 0;
            comsolA = 0;
            nsolA = 0;
        else
        R2 = [1 0 0; 0 cos(-t1s2) -sin(-t1s2); 0 sin(-t1s2) cos(-t1s2)];
        center1 = R2*[0;r;0];
        midpoint_dist = midpoints - repmat(center1,[1,size(midpoints,2)]);
        midpoint_dist = midpoint_dist.'*midpoint_dist;
        midpoint_dist = diag(midpoint_dist);
        ind = find(midpoint_dist==min(midpoint_dist));
        midpoint1 = midpoints(:,ind);
        [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC2(N,r,p1,p2,p3,t1s2,midpoint1);
        end
    case 0
        disp('0');
        midpoints = TS_sidefingeredge(r,p1,p2,p3,phi);
        
        
        if isempty(midpoints)
            cliffHeight = inf;
            centersolA = 0;
            comsolA = 0;
            nsolA = 0;
        else
        end
end
end
% toc
if draw
%% draw figures for paper
% iterated cliff curve
specified = 350;
[fig,h] = draw3DSpheres(r,p1,p2,p3);

h = drawIterationResult(fig,h,r,p1,p2,p3,centersolA,comsolA,specified,iterate,'k');
result = [centersolA(specified,:).' comsolA(specified,:).'];
h(end+1) = plot3(result(1,:),result(2,:),result(3,:),'-+k','lineWidth',3,'displayName','iter. 1');
iterate = 2;
[cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC3(phi3_rangeA,r,p1,p2,p3,t1s1,t1s2,iterate);
h = drawIterationResult(fig,h,r,p1,p2,p3,centersolA,comsolA,specified,iterate,'--g');
result = [centersolA(specified,:).' comsolA(specified,:).'];
h(end+1) = plot3(result(1,:),result(2,:),result(3,:),'-+g','lineWidth',3,'displayName','iter. 2');
iterate = 10;
tic
[cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC3(phi3_rangeA,r,p1,p2,p3,t1s1,t1s2,iterate);
toc
h = drawIterationResult(fig,h,r,p1,p2,p3,centersolA,comsolA,specified,iterate,'y');
result = centersolA([1,end/2],:).';
h(end+1) = plot3(result(1,:),result(2,:),result(3,:),'+r','lineWidth',4,'markerSize',4,'displayName','symm. points');
legend(h);


%cone process
specified = 70;
[fig,h] = draw3DSpheres(r,p1,p2,p3);
iterate = 1;
[cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC3(phi3_rangeA,r,p1,p2,p3,t1s1,t1s2,iterate);
drawIterationResult(fig,h,r,p1,p2,p3,centersolA,comsolA,specified,iterate,'k');
h = drawCones(fig,h,r,p1,centersolA(specified,:));
h = drawCones(fig,h,r,p2,centersolA(specified,:));
result = [centersolA(specified,:).' comsolA(specified,:).'];
h(end) = plot3(result(1,:),result(2,:),result(3,:),'k','lineWidth',4,'displayName','cyl. axis');
h(end+1) = plot3(result(1,1),result(2,1),result(3,1),'+g','lineWidth',2,'markerSize',4,'displayName','base center');
h(end+1) = plot3(result(1,2),result(2,2),result(3,2),'+r','lineWidth',2,'markerSize',6,'displayName','COM');
legend(h);
end
switch nargout
    case 0
        % save('cliff4.mat','phi3_range','cliff_height','cliffTheta1','cliffTheta2');
        figure
        hold on
        plot(phi3_range,cliffHeight)
        % plot(cliffTheta1,cliff_height)
        xlabel('\phi_3')
        ylabel('COM height');
%         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,minCOM,minTheta1,minTheta2,minPhi3,p1,p2,p3);
    fig = draw3DSpheres(r,p1,p2,p3);
    hold on
    
    for i=1:size(nsolA,1)
        centerA = centersolA(i,:).';
        comA = comsolA(i,:).';
%         axisA = [centerA comA];
        axisA = [centerA];
        plot3(axisA(1,:),axisA(2,:),axisA(3,:),'.k','lineWidth',4,'markerSize',2);
    end
    hold off
    
    case 1
        varargout{1} = cliffHeight;
    case 2
        varargout{1} = phi3_rangeA;
        varargout{2} = cliffHeight;
    case 4
        varargout{1} = phi3_rangeA;
        varargout{2} = cliffHeight;
        varargout{3} = nsolA;
        varargout{4} = comsolA;
    case 5
        varargout{1} = phi3_rangeA;
        varargout{2} = cliffHeight;
        varargout{3} = nsolA;
        varargout{4} = comsolA;
        varargout{5} = centersolA;
%     case 3
%         varargout{1} = minCOM;
%         varargout{2} = minTheta1;
%         varargout{3} = minTheta2;
end
end

function h = drawIterationResult(fig,h,r,p1,p2,p3,centersolA,comsolA,specified,iter,col)
set(0, 'currentfigure', fig);
hold on
    
    result = centersolA.';
    h(end+1) = plot3(result(1,:),result(2,:),result(3,:),col,'lineWidth',2,'displayName',['iter. ' num2str(iter)]);
    if isempty(specified)
        specified = 0;
        
        for i=1:size(centersolA,1)
            comA = comsolA(i,:).';
            %         axisA = [centerA comA];
            
            
            if i==specified
                result = [centerA comA];
                plot3(result(1,:),result(2,:),result(3,:),'k','lineWidth',2,'markerSize',2);
            end
        end
    end
end
function h = drawCones(fig,h,r,p1,centersol)
set(0, 'currentfigure', fig);
d1 = p1-centersol.';
nd1 = norm(d1);
% d1n = d1/nd1;
ph1 = asin(r/norm(d1));
hold on
%compute angle between d1 and X-Y plane
alpha1 = asin(d1(3)/nd1);
%compute angle revolving d1 around z-axis, relative to x-axis.
alpha2 = atan2(d1(2),d1(1));
[X,Y,Z]=cylinder([0 tan(ph1)*nd1]);
Z = Z*nd1;
Caxis = [centersol.' centersol.'+d1];
% plot3(Caxis(1,:),Caxis(2,:),Caxis(3,:),'k','lineWidth',4,'markerSize',2);
M=makehgtform('translate',centersol,'zrotate',alpha2,'yrotate',pi()/2-alpha1);
h(end+1)=surf(X,Y,Z,'Parent',hgtransform('Matrix',M),'FaceAlpha',0.4,'FaceColor','g','displayName','cones');%,'LineStyle','none'
view([30,35])
grid on
light       
        
end

function [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC3(phi3_range,r,p1,p2,p3,t1s1,t1s2,iterate)
cliffHeight = inf*ones(size(phi3_range));
centersolA = inf*ones(size(phi3_range,1),3);
comsolA = centersolA;
nsolA = centersolA;
p2t = p2-p3;
p1t = p1-p3;


R1 = [1 0 0; 0 cos(-t1s1) -sin(-t1s1); 0 sin(-t1s1) cos(-t1s1)];
R2 = [1 0 0; 0 cos(-t1s2) -sin(-t1s2); 0 sin(-t1s2) cos(-t1s2)];
t1c1center = R1*[0;-r;0]; %calculated relative to p3
t1c2center = R2*[0;r;0];
circleCenter = (t1c1center+t1c2center)/2;
noon = t1c1center-circleCenter; %a radius vector for the circle which serves as the initial guess for the center curve

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
    for iter = 1:iterate
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
        if iter<iterate
            center = pointOnLine*r/norm(pointOnLine);
        end
    end
    
    %some plots to look at certain values relevant to the correction
%     plot(i,r-norm(pointOnLine),'g+')
%     plot(i,sol*center,'r+')
    


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

function [cliffHeight,centersolA,comsolA,nsolA] = computeHalfCliffApprox2ndC2(N,r,p1,p2,p3,t1s2,center2)
cliffHeight = inf*ones(N,1);
centersolA = inf*ones(N,3);
comsolA = centersolA;
nsolA = centersolA;
p2t = p2-p3;
p1t = p1-p3;
%% find other endpoint of the cliff curve section:
% in this case it is the point in which O1 reaches the base edge:

%%
R2 = [1 0 0; 0 cos(-t1s2) -sin(-t1s2); 0 sin(-t1s2) cos(-t1s2)];
center1 = R2*[0;r;0];
angle_range = acos(center1.'*center2/(norm(center1)*norm(center2)));
normal = cross(center1,center2);
u = cross(normal,center1); %u is vector normal to center1 in the plane spanned by center1,center2
u = r*u/norm(u);
fig = draw3DSpheres(r,p1,p2,p3);
hold on
plot3(center1(1),center1(2),center1(3),'+g');
plot3(center2(1),center2(2),center2(3),'+g');

for i=1:N
    phi = (i-1)/(N-1)*angle_range;
    center = center1*cos(phi)+u*sin(phi);
    
    % check if the base center is within the sphere of one of the other
    % fingers: if yes, move the calculation
    d1 = p1t-center;
    nd1 = norm(d1);
    d2 = p2t-center;
    nd2 = norm(d2);
    if nd1<r
        normal = cross(p1-p3,center);
        normal = normal/norm(normal);
        p1tn = p1t/norm(p1t);
        utag = cross(p1tn,normal);
        utag = utag/norm(utag); %this vector is perpendicular to the line connecting p1,p3, and lies in the same plane as center
        center_rel = (center-0.5*(p1+p3));
        ksi = acos(center_rel.'*utag/norm(center_rel));
        center = 0.5*(p1+p3)+norm(center_rel)*(utag*cos(ksi)-p1tn*sin(ksi));
    end
    centerp = center+p3;
    plot3(centerp(1),centerp(2),centerp(3),'+k');
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
    plot3(center(1),center(2),center(3),'+r');
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



function base_center = TS_sidefingeredge(r,p1,p2,p3,phi)
%function returns the base center relative to p3, for the stance in which
%O3,O1 touch base edge
base_center = [];
%the intersection of the spheres surrounding p1,p3 is the following circle:
circle1_center = 0.5*(p1+p3);
circle1_radius = sqrt(r^2-0.25*norm(p1-p3)^2);
%the plane defined by this circle intersects the sphere surrounding p2 in
%the following circle:
N = p1-p3; %normal to the plane
N = N/norm(N);
D = -circle1_center.'*N;
dist_p2_to_plane = abs(N.'*p2+D);
if dist_p2_to_plane>r
    return
end
circle2_center = p2+dist_p2_to_plane*N;
circle2_radius = sqrt(r^2-dist_p2_to_plane^2);
%% the lines tangent to both circles are as following:

relcc = circle2_center-circle1_center;
relcc_norm = norm(relcc);
fig = draw3DSpheres(r,p1,p2,p3);
hold on
plot3(circle1_center(1),circle1_center(2),circle1_center(3),'r+','markersize',8,'linewidth',6);
plot3(circle2_center(1),circle2_center(2),circle2_center(3),'r+','markersize',8,'linewidth',6);
phi_range = 0:0.01:2*pi();
for i=1:length(phi_range)
    phi = phi_range(i);
    u = relcc/norm(relcc);
    v = cross(u,N);
    circle1(:,i) = circle1_center+circle1_radius*(u*cos(phi)+v*sin(phi));
end
plot3(circle1(1,:),circle1(2,:),circle1(3,:),'r','linewidth',4);
for i=1:length(phi_range)
    phi = phi_range(i);
    u = relcc/norm(relcc);
    v = cross(u,N);
    circle2(:,i) = circle2_center+circle2_radius*(u*cos(phi)+v*sin(phi));
end
plot3(circle2(1,:),circle2(2,:),circle2(3,:),'r','linewidth',4);
if relcc_norm<circle1_radius+circle2_radius %if the two circles intersect:
    % there are two possible solutions
    if circle2_radius>circle1_radius
        x = relcc_norm*circle1_radius/(circle2_radius-circle1_radius);
        angle_axis_relcc = asin(circle1_radius/x);
        angle_base_center = pi()/2-angle_axis_relcc;
        % construct the radius tangent to the axis from circle1_center:
        u = -relcc/relcc_norm;
        v = cross(N,u);
        base_center(:,1) = circle1_center -p3 + circle1_radius*(u*cos(angle_base_center)+v*sin(angle_base_center));
        base_center(:,2) = circle1_center -p3 + circle1_radius*(u*cos(angle_base_center)-v*sin(angle_base_center));
    else
        x = relcc_norm*circle2_radius/(circle1_radius-circle2_radius);
        angle_axis_relcc = asin(circle2_radius/x);
        angle_base_center = pi()/2-angle_axis_relcc;
        % construct the radius tangent to the axis from circle1_center:
        u = relcc/relcc_norm;
        v = cross(N,u);
        base_center(:,1) = circle1_center -p3 + circle1_radius*(u*cos(angle_base_center)+v*sin(angle_base_center));
        base_center(:,2) = circle1_center -p3 + circle1_radius*(u*cos(angle_base_center)-v*sin(angle_base_center));
    end
    norm(base_center(:,1))
    %     axis = relcc + base_center*(circle2_radius/circle1_radius-1);
    %     if axis(3)<0
    %         disp('chose wrong solution in 2 solution case');
    %     end
else %the two circles do not intersect
    %there are 4 possible solutions
    
end

end
