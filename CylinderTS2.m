function [cliff_phi3,cliff_s3,cliff_height] = CylinderTS2(r,hcm,P,t,phi,Ni,print)
s3range = linspace(0,r,Ni);%linspace(0,r,N);%
% s3range = [s3range linspace(0.71*r,r,2*Ni)];
N = length(s3range);
Nphi = 50;
heights = -inf*ones(N,2*Nphi);
[t_start,~] = TS_symmetrical2(r,0,P,phi);
prev_sol = [t_start*ones(2,Nphi)];
%%

[feq3,t2_sol3] = generate_feq_fast(r,P);


%%
countinf = 0;
% phi3_range = linspace(-pi(),0,Nphi);
for i=1:N
    s3 = s3range(i)
    %     feq2 = @(x,phi3) feq3(x,phi3,s3);
    t1_feq31 = feq3{1};
    t1_feq32 = feq3{2};
    t1_feq21 = @(x,phi3) t1_feq31(x,phi3,s3);
    t1_feq22 = @(x,phi3) t1_feq32(x,phi3,s3);
    feq2 = {t1_feq21,t1_feq22};
    t2_sol2 = @(x,phi3) t2_sol3(x,phi3,s3);
%     tic
    [~,height,theta1,theta2] = CylinderTSCurve3(r,s3,hcm,P,t,phi,Nphi,prev_sol,feq2,t2_sol2);
%     toc_it = toc
    prev_sol = [theta1;theta2];
    height(height==inf)=-inf;
    if all(height==-inf)
        countinf = countinf + 1;
        if countinf>3
            break
        end
    end
    heights(i,1:Nphi)=height;
    heights(i,Nphi+1:2*Nphi)=height(Nphi:-1:1);

end
phi3_range = linspace(-pi(),0,Nphi);
cliff_phi3 = phi3_range;
cliff_s3 = phi3_range;
cliff_height = phi3_range;
for j=1:length(phi3_range)
    height_ind = find(heights(:,j)==-inf,1,'first')-1;
    if isempty(height_ind)
        height_ind = length(s3range);
    end
    cliff_s3(j) = s3range(height_ind);
    cliff_height(j) = heights(height_ind,j);
end

if print>0
    threeD = print==3;
fig1 = figure;
cliff_height = [cliff_height cliff_height(Nphi:-1:1)];
cliff_s3 = [cliff_s3 cliff_s3(Nphi:-1:1)];
cliff_phi3 = [cliff_phi3 -cliff_phi3(Nphi:-1:1)];
phi3_range = [phi3_range,-phi3_range(Nphi:-1:1)];
try
    [s3range_grid,phi3_range_grid] = meshgrid(s3range,phi3_range);
    x = s3range_grid.*sin(phi3_range_grid);
    y = s3range_grid.*cos(phi3_range_grid);
%     x = [x,-x];
%     y = [y,y];
%     heights = [heights;heights];
    % [~,cf] = contourf(s3range,phi3_range,heights.',30);
    if threeD
        sur = surf(x,y,heights.');
        zlabel('$z$','interpreter','latex','FontSize',20)
        set(get(gca,'zlabel'),'rotation',0);
        set(sur,'EdgeAlpha',0.1)
    else
     [~,cf] = contourf(x,y,heights.',30);
     cb = colorbar;
     set(get(cb,'label'),'string','$\:U_{3c}$','interpreter','latex','FontSize',20,'rotation',0);
    end
    
    hold on

    % 
axis equal
% xlabel('$s_3$','interpreter','latex','FontSize',20)
% ylabel('$\phi_3$','interpreter','latex','FontSize',20)
xlabel('$\tilde{x}$','interpreter','latex','FontSize',20)
ylabel('$\tilde{y}$','interpreter','latex','FontSize',20)

%     plot3(cliff_x,cliff_y,cliff_height,'r','lineWidth',6)
catch me
    disp('oops');
end



end
end

function feq = generate_feq(r,P)
x = sym('x',[2,1]);
t1t = x(1);
t2t = x(2);
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
syms phi3 d3
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rt*Rz*[0 -d3 0].';
n = Rt*nt;
p1p = p1t-(p1t.'*n)*n;
p2p = p2t-(p2t.'*n)*n;
rad1 = p1p-center;
rad2 = p2p-center;
eq1 = rad1.'*rad1-r^2;
eq2 = rad2.'*rad2-r^2;
feq = matlabFunction([eq1;eq2],'vars',{x,phi3,d3});
end
function [feq,t2_feq] = generate_feq2(r,P)
syms t2t
x = sym('x',[2,1]);
t1t = x(1);
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
syms phi3 d3
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rt*Rz*[0 -r 0].';
n = Rt*nt;
rad1 = cross(p1t-center,n);
eq1 = rad1.'*rad1-r^2;
% eq2 = rad2.'*rad2-r^2;
t2_sol = [-asin((r*sin(phi3))/(p1(3)*cos(x(1)) - p3(3)*cos(x(1)) - p1(2)*sin(x(1)) + p3(2)*sin(x(1))));...
 pi + asin((r*sin(phi3))/(p1(3)*cos(x(2)) - p3(3)*cos(x(2)) - p1(2)*sin(x(2)) + p3(2)*sin(x(2))))];
eq(1,1) = subs(eq1,t2t,t2_sol(1));
eq(2,1) = subs(eq1,[t2t,x(1)],[t2_sol(2),x(2)]);
feq = matlabFunction(eq,'vars',{x,phi3,d3});
t2_feq = matlabFunction(t2_sol,'vars',{x,phi3,d3});
end
function [feq,t2_feq] = generate_feq3(r,P)
syms t2t
x = sym('x',[2,1]);
t1t = x(1);
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p1y = p1(2);
p1z = p1(3);
p3y = p3(2);
p3z = p3(3);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
syms phi3 d3 t1
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rt*Rz*[0 -d3 0].';
n = Rt*nt;
rad1 = cross(p1t-center,n);
eq1 = rad1.'*rad1-r^2;
% eq2 = rad2.'*rad2-r^2;
t2_sol = [-asin((d3*sin(phi3))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))));...
 pi + asin((d3*sin(phi3))/(p1z*cos(x(2)) - p3z*cos(x(2)) - p1y*sin(x(2)) + p3y*sin(x(2))))];
eq(1,1) = subs(eq1,t2t,t2_sol(1))
eq(2,1) = subs(eq1,[t2t,x(1)],[t2_sol(2),x(2)])
eq = subs(eq,[x(1),x(2)],[t1,t1]);
% feq = matlabFunction(eq,'vars',{x,phi3,d3});
t1_feq1 = matlabFunction(eq(1,1),'vars',{t1,phi3,d3});
t1_feq2 = matlabFunction(eq(2,1),'vars',{t1,phi3,d3});
feq = {t1_feq1,t1_feq2};
% feq = @(x,phi,d3) [
t2_feq = matlabFunction(t2_sol,'vars',{x,phi3,d3});
end
function [t1_feq,t2_feq] = generate_feq_fast(r,P)
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p1x = p1(1);
p1y = p1(2);
p1z = p1(3);
p3y = p3(2);
p3z = p3(3);
t2_feq = @(x,phi3,d3) -asin((d3*sin(phi3))./(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))));
% t1_feq = {@(x,phi3,d3) (cos(x(1))*(p1x - d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) + (d3*sin(phi3)*(p1z - p3z + d3*cos(phi3)*sin(x(1)) - (d3^2*cos(x(1))*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2 - r^2 - ((d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2 - 1)*(d3*cos(phi3) + p1y*cos(x(1)) - p3y*cos(x(1)) + p1z*sin(x(1)) - p3z*sin(x(1)))^2 + (sin(x(1))*(p1x - d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) - (d3*sin(phi3)*(p1y - p3y + d3*cos(phi3)*cos(x(1)) + (d3^2*sin(phi3)^2*sin(x(1)))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2,...
%  @(x,phi3,d3) (cos(x(1))*(p1x + d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) - (d3*sin(phi3)*(p1z - p3z + d3*cos(phi3)*sin(x(1)) - (d3^2*cos(x(1))*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2 - r^2 - ((d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2 - 1)*(d3*cos(phi3) + p1y*cos(x(1)) - p3y*cos(x(1)) + p1z*sin(x(1)) - p3z*sin(x(1)))^2 + (sin(x(1))*(p1x + d3*sin(phi3)*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2))*(1 - (d3^2*sin(phi3)^2)/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))^2)^(1/2) + (d3*sin(phi3)*(p1y - p3y + d3*cos(phi3)*cos(x(1)) + (d3^2*sin(phi3)^2*sin(x(1)))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1)))))/(p1z*cos(x(1)) - p3z*cos(x(1)) - p1y*sin(x(1)) + p3y*sin(x(1))))^2};
t1_feq = {@(x,phi3,d3) (cos(x).*(p1x - d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) + (d3.*sin(phi3).*(p1z - p3z + d3.*cos(phi3).*sin(x) - (d3.^2.*cos(x).*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2 - r.^2 - ((d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2 - 1).*(d3.*cos(phi3) + p1y.*cos(x) - p3y.*cos(x) + p1z.*sin(x) - p3z.*sin(x)).^2 + (sin(x).*(p1x - d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) - (d3.*sin(phi3).*(p1y - p3y + d3.*cos(phi3).*cos(x) + (d3.^2.*sin(phi3).^2.*sin(x))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2,...
 @(x,phi3,d3) (cos(x).*(p1x + d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) - (d3.*sin(phi3).*(p1z - p3z + d3.*cos(phi3).*sin(x) - (d3.^2.*cos(x).*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2 - r.^2 - ((d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2 - 1).*(d3.*cos(phi3) + p1y.*cos(x) - p3y.*cos(x) + p1z.*sin(x) - p3z.*sin(x)).^2 + (sin(x).*(p1x + d3.*sin(phi3).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2)).*(1 - (d3.^2.*sin(phi3).^2)./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x)).^2).^(1./2) + (d3.*sin(phi3).*(p1y - p3y + d3.*cos(phi3).*cos(x) + (d3.^2.*sin(phi3).^2.*sin(x))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))))./(p1z.*cos(x) - p3z.*cos(x) - p1y.*sin(x) + p3y.*sin(x))).^2};
end
