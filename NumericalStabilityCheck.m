function [flag,val] = NumericalStabilityCheck(r,p1,p2,p3,t,N)

phi3_range = linspace(-pi(),0,N*20);
[CurveHeight,centersolA,comsolA,nsolA]  = computeCurveAroundGrasp(phi3_range,r,p1,p2,p3,t); %#ok<ASGLU>

minHeight = min(CurveHeight);
% GraspHeight = CylinderComPos(r,1,p1,p2,p3,t,0,s3);
GraspHeight = 0;%cos(t);
LocalDepth = minHeight - GraspHeight;
flag = (sign(LocalDepth)+1)/2;
val = LocalDepth;

% figure
% hold on
% plot3(centersolA(:,1),centersolA(:,2),centersolA(:,3))
% plot3(0,0,0,'g+')
% figure
% hold on
% plot3(comsolA(:,1),comsolA(:,2),comsolA(:,3))
% plot3(0,-sin(t),cos(t),'g+')
% f = figure;
% hold on
% plot3(p1(1),p1(2),p1(3),'or','markerSize',4,'lineWidth',4);
% plot3(p2(1),p2(2),p2(3),'ob','markerSize',4,'lineWidth',4);
% plot3(p3(1),p3(2),p3(3),'og','markerSize',4,'lineWidth',4);
% [x,y,z] = sphere;
% mesh(r*x+p1(1),r*y+p1(2),r*z+p1(3),'faceAlpha',0.2);
% mesh(r*x+p2(1),r*y+p2(2),r*z+p2(3),'faceAlpha',0.2);
% mesh(r*x+p3(1),r*y+p3(2),r*z+p3(3),'faceAlpha',0.2);
% 
% for i=1:size(nsolA,1)
%     
%     centerA = centersolA(i,:).';
%     comA = comsolA(i,:).';
%     axisA = [centerA comA];
%     
%     plot3(axisA(1,:),axisA(2,:),axisA(3,:),'ro-','lineWidth',4,'markerSize',2);
% end
% axis equal
% xlabel('x')
% ylabel('y')
% zlabel('z')
end
% function minHeight = computeCurveAroundGrasp(phi3_range,r,p1,p2,p3,t)

function [CurveHeight,centersolA,comsolA,nsolA] = computeCurveAroundGrasp(phi3_range,r,p1,p2,p3,t)

p2t = p2;
p1t = p1;


noon = [0;cos(t);sin(t)];
% circleAxis = [0;-sin(tmean);cos(tmean)];
% figure
% hold on
perpToBoth = [-1;0;0];
% scalev = linspace(0,r/10,100);%logspace(-5,log(r/10)/log(10),50);
scalev = 0.01;
for j=1:length(scalev)
    scale = scalev(j);
    CurveHeight = inf*ones(size(phi3_range));
centersolA = inf*ones(size(phi3_range,2),3);
comsolA = centersolA;
nsolA = centersolA;
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
%     Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
    clockArm = scale*(noon*cos(phi3) + perpToBoth*sin(phi3));
    centerb = clockArm; %rotating the centerbase around the origin, which is where the center base is at the BG

    %compute parameters of cones starting from center tangent to spheres
    %surrounding O1,2
    d1 = p1t-centerb;
    nd1 = norm(d1);
    d1n = d1/nd1;
    ph1 = asin(r/nd1);
    
    d2 = p2t-centerb;
    nd2 = norm(d2);
    d2n = d2/nd2;
    ph2 = asin(r/nd2);
    %compute COM position given tangency criterion+base center position
    coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
    caxis = [0;-sin(t);cos(t)]; % axis at the basket grasp
    if j==1 && coms(1,:)*caxis>coms(2,:)*caxis
        %take the solution closer to the direction of the axis at the BG
        sol_choice = 1;
    else
        sol_choice = 2;
    end
    sol = coms(sol_choice,:);
%     if ~imag(sol)==0
% %         continue;
%         if imag(coms(1,:))==0
%         sol = coms(1,:);
%         else
%             continue;
%         end
%     end
    pointOnLine1 = -(centerb-p3).'*sol.'/(sol*sol.')*sol.'+centerb;
    
%     sol*center
%     norm(pointOnLine)
    centerb = pointOnLine1;
%     dist = (centerb-p3).'*sol.'
    %some plots to look at certain values relevant to the correction
%     plot(i,r-norm(pointOnLine),'g+')
%     plot(i,sol*center,'r+')
    
%%additional correction  - rerun the computation of the tangent using
    %%the corrected center - turned out to get worse results, sus?
%     d1 = p1t-centerb;
%     d1n = d1/norm(d1);
%     ph1 = asin(r/norm(d1));
%     
%     d2 = p2t-centerb;
%     d2n = d2/norm(d2);
%     ph2 = asin(r/norm(d2));
%     
%     coms = comGivenBaseCenter(d1n,d2n,ph1,ph2);
%     sol = coms(2,:); 
%     
%     pointOnLine = -(centerb-p3).'*sol.'/(sol*sol.')*sol.'+centerb;
%     centerb = pointOnLine;
% 
%     change = norm(pointOnLine-pointOnLine1)


    com = sol.'+centerb;
    nsol = sol.'/norm(sol);
    
    nsolA(i,:) = nsol.';
    centersolA(i,:) = centerb.';
    comsolA(i,:) = com.';
    ht = com(3);
    CurveHeight(i) = ht-cos(t); 
    
    
end
heights(j,:) = CurveHeight;


end
% minHeight = min(CurveHeight);
% fig1 = figure;
% 
%     try
%     [s3range_grid,phi3_range_grid] = meshgrid(scalev,phi3_range);
%     x = s3range_grid.*sin(phi3_range_grid);
%     y = s3range_grid.*cos(phi3_range_grid);
% %     x = [x,-x];
% %     y = [y,y];
% %     heights = [heights;heights];
%     % [~,cf] = contourf(s3range,phi3_range,heights.',30);
%     [~,cf] = contourf(x,y,heights.',30);
% catch me
% %     disp('oops');
% end
% cb = colorbar;
% axis equal
% xlabel('$x$','interpreter','latex','FontSize',20)
% ylabel('$y$','interpreter','latex','FontSize',20)
% set(get(cb,'label'),'string','$\:U_{3c}$','interpreter','latex','FontSize',20,'rotation',0);
end
