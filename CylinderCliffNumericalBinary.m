function [cliff_phi3,cliff_s3,cliff_height,comsolA,nsolA] = CylinderCliffNumericalBinary(r,hcm,t,P,Ni,tolerance,print)
%linspace(0,r,N);%
% s3range = [s3range linspace(0.71*r,r,2*Ni)];

options0 = optimset('Display','none','TolX',1e-8);

%%
t1_feq3 = generate_feq_fast(r,P);


%%
s3 = norm(P(:,3));

cliff_phi3 = linspace(-pi(),0,Ni);
cliff_s3 = cliff_phi3;

cliff_height = cliff_phi3;
comsolA = zeros(Ni,3);
anglesA = zeros(Ni,2);
nsolA = zeros(Ni,3);
for i=1:Ni
    %first a computation to find the upper bound for the distance of the
    %cliff curve point at this angle from the basket grasp.
    found_height = inf;
    phi3_rel = cliff_phi3(i);
    a = tan(phi3_rel);
    ys = roots([1+a^2,2*s3,s3^2-r^2]);
%     y = ys(abs(ys)<=r);
    signy = sign(r*cos(phi3_rel));
    if length(ys)>1
        yrel = ys(sign(ys) == signy);
    else
        disp('single solution in range');
    end
    if isempty(yrel)
        disp('no solution in range');
    end
    y = yrel+s3;
    if abs(y)>r
        num_err = abs(y)-r
        y = sign(y)*r;
    end
    x = -sqrt(r^2-y^2);
    direct = [x;y]-[0;s3];
    if abs(wrapToPi(atan2(direct(1),direct(2)))-phi3_rel)>1E-6
        disp('CCNB: error in boundary computation')
    end
    directn = direct/norm(direct);
    bounds = [0,norm(direct)];
    
    
    found_cliff = 0;
    while ~found_cliff
%     tic
    d3_rel = mean(bounds);
    p3_obj = [0;s3]+d3_rel*directn;
    d3 = norm(p3_obj);
    phi3 = atan2(p3_obj(1),p3_obj(2));
    t1_feq1 = @(x) t1_feq3(x,phi3,d3);
    [height,comsol,angles,nsol] = Triple_support_stance_numerical(r,hcm,t,P,t1_feq1,phi3,d3,options0); %#ok<ASGLU>
%     [~,height,theta1,theta2] = CylinderTSCurve3(r,s3,hcm,P,t,phi,Nphi,prev_sol,feq2,t2_sol2);
%     toc_it = toc %#ok<NOPRT,NASGU>
    if height==inf
        bounds = [bounds(1),d3_rel];
    else
        prev_height = found_height;
        found_height = height;
        found_comsol = comsol;
        found_angles = angles;
        found_nsol = nsol;
        found_phi3 = phi3;
        found_d3 = d3;
        bounds = [d3_rel,bounds(2)];
        if abs(prev_height-found_height)<tolerance
            found_cliff = 1;
        end
    end
%     if diff(bounds)<tolerance/100&&found_height==inf
%         bounds = [-norm(P(:,3)),0];
%     end
    end
%     if found_height == inf
%         [found_height,found_comsol] = Triple_support_stance_numerical(r,hcm,t,P,t1_feq2,phi3,d3,options0);
%         disp('no stance found for specific phi3');
%         bounds = [0,0];
%        
%     end
    cliff_height(i) = found_height;
    cliff_s3(i) = found_d3;
    cliff_phi3(i) = found_phi3;
    comsolA(i,:) = found_comsol.';
    anglesA(i,:) = found_angles;
    nsolA(i,:) = found_nsol.';
%     if all(height==-inf)
%         countinf = countinf + 1;
%         if countinf>3
%             break
%         end
%     end


end



if print

cliff_s3_p = [cliff_s3 cliff_s3(Ni:-1:1)];
cliff_phi3_p = [cliff_phi3 -cliff_phi3(Ni:-1:1)];

cliff_x = cliff_s3_p.*sin(cliff_phi3_p);
cliff_y = cliff_s3_p.*cos(cliff_phi3_p);
hold on
plot(cliff_x,cliff_y,'k','lineWidth',4)


hold off
axis equal
% % xlabel('$s_3$','interpreter','latex','FontSize',20)
% % ylabel('$\phi_3$','interpreter','latex','FontSize',20)
% xlabel('$x$','interpreter','latex','FontSize',20)
% ylabel('$y$','interpreter','latex','FontSize',20)
end
end


function [t1_feq] = generate_feq_fast(r,P)
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p1x = p1(1);
p1y = p1(2);
p1z = p1(3);
p3y = p3(2);
p3z = p3(3);
% t2_feq = @(x,phi3,d3) -asin((d3*sin(phi3))./(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))));
% t1_feq = {@(x,phi3,d3) (cos(x(1))*(p1x - d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) + (d3*sin(phi3)*(p1z - p3z + d3*cos(phi3)*sin(x(1)) - (d3^2*cos(x(1))*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2 - r^2 - ((d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2 - 1)*(d3*cos(phi3) + p1y*cos(x(1)) - p3y*cos(x(1)) + p1z*sin(x(1)) - p3z*sin(x(1)))^2 + (sin(x(1))*(p1x - d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) - (d3*sin(phi3)*(p1y - p3y + d3*cos(phi3)*cos(x(1)) + (d3^2*sin(phi3)^2*sin(x(1)))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2,...
%  @(x,phi3,d3) (cos(x(1))*(p1x + d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) - (d3*sin(phi3)*(p1z - p3z + d3*cos(phi3)*sin(x(1)) - (d3^2*cos(x(1))*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2 - r^2 - ((d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2 - 1)*(d3*cos(phi3) + p1y*cos(x(1)) - p3y*cos(x(1)) + p1z*sin(x(1)) - p3z*sin(x(1)))^2 + (sin(x(1))*(p1x + d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) + (d3*sin(phi3)*(p1y - p3y + d3*cos(phi3)*cos(x(1)) + (d3^2*sin(phi3)^2*sin(x(1)))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2};
t1_feq = @(x,phi3,d3) (cos(x).*(p1x - d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) + (d3.*sin(phi3).*(p1z - p3z + d3.*cos(phi3).*sin(x) - (d3.^2.*cos(x).*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2 - r.^2 - ((d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2 - 1).*(d3.*cos(phi3) + p1y.*cos(x) - p3y.*cos(x) + p1z.*sin(x) - p3z.*sin(x)).^2 + (sin(x).*(p1x - d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) - (d3.*sin(phi3).*(p1y - p3y + d3.*cos(phi3).*cos(x) + (d3.^2.*sin(phi3).^2.*sin(x))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2;
%  @(x,phi3,d3) (cos(x).*(p1x + d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) - (d3.*sin(phi3).*(p1z - p3z + d3.*cos(phi3).*sin(x) - (d3.^2.*cos(x).*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2 - r.^2 - ((d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2 - 1).*(d3.*cos(phi3) + p1y.*cos(x) - p3y.*cos(x) + p1z.*sin(x) - p3z.*sin(x)).^2 + (sin(x).*(p1x + d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) + (d3.*sin(phi3).*(p1y - p3y + d3.*cos(phi3).*cos(x) + (d3.^2.*sin(phi3).^2.*sin(x))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2};
end
