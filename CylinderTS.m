function CylinderTS(r,hcm,P,t,phi,Ni)
s3range = linspace(0,r,Ni);%linspace(0,r,N);%
% s3range = [s3range linspace(0.71*r,r,2*Ni)];
N = length(s3range);
Nphi = 150;
heights = -inf*ones(N,2*Nphi);
[t_start,~] = TS_symmetrical2(r,0,P,phi);
prev_sol = [t_start*ones(1,Nphi);zeros(1,Nphi)];
%%
[feq3,t2_sol] = generate_feq(r,P);
% t2_sol = subs(t2_sol)
%%
countinf = 0;
% phi3_range = linspace(-pi(),0,Nphi);
for i=1:N
    s3 = s3range(i)
    feq2 = @(x,phi3) feq3(x,phi3,s3);
    [phi3_range,height,theta1,theta2] = CylinderTSCurve(r,s3,hcm,P,t,phi,Nphi,prev_sol,feq2,t2_sol);
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
for j=1:length(phi3_range)
    cliff_phi3(j) = phi3_range(j);
    height_ind = find(heights(:,j)==-inf,1,'first')-1;
    cliff_s3(j) = s3range(height_ind);
    cliff_height(j) = heights(height_ind,j);
end


fig1 = figure;
phi3_range = [phi3_range,-phi3_range(Nphi:-1:1)];
try
    [s3range_grid,phi3_range_grid] = meshgrid(s3range,phi3_range);
    x = s3range_grid.*sin(phi3_range_grid);
    y = s3range_grid.*cos(phi3_range_grid);
%     x = [x,-x];
%     y = [y,y];
%     heights = [heights;heights];
    % [~,cf] = contourf(s3range,phi3_range,heights.',30);
    [~,cf] = contourf(x,y,heights.',30);
    hold on
    cliff_x = cliff_s3.*sin(cliff_phi3);
    cliff_y = cliff_s3.*cos(cliff_phi3);
    plot(cliff_x,cliff_y,'k')
catch me
    disp('oops');
end


cb = colorbar;
axis equal
% xlabel('$s_3$','interpreter','latex','FontSize',20)
% ylabel('$\phi_3$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
ylabel('$y$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U_{3c}$','interpreter','latex','FontSize',20,'rotation',0);
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
function [feq,t2_sol] = generate_feq2(r,P)
syms x t2t
t1t = x;
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p2t = p2-p3;
p1t = p1-p3;
Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
Rt = Rx*Ry;
nt = [0;0;1];
syms phi3
Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
center = Rt*Rz*[0 -d3 0].';
n = Rt*nt;
p1p = p1t-(p1t.'*n)*n;
p2p = p2t-(p2t.'*n)*n;
rad1 = cross(p1t-center,n);
eq = rad1.'*rad1-r^2;
% eq2 = rad2.'*rad2-r^2;
t2_sol = [-asin((r*sin(phi3))/(p1(3)*cos(x) - p3(3)*cos(x) - p1(2)*sin(x) + p3(2)*sin(x)));...
 pi + asin((r*sin(phi3))/(p1(3)*cos(x) - p3(3)*cos(x) - p1(2)*sin(x) + p3(2)*sin(x)))];
eq(1,1) = subs(eq,t2t,t2_sol(1));
eq(2,1) = subs(eq,t2t,t2_sol(2));
feq = matlabFunction(eq,'vars',{x1,phi3,d3});
end
