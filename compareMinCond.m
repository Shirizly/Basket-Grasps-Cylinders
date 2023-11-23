% syms t p1 p2
for i = 1
t = (rand(1))*pi()/2;
ph = pi()/2*(rand(1)+1);
ph1 = ph; ph2 = -ph;
r = 0.5;
h = 2;
hcm = 0.5*h;
com = [0;0;hcm];
ds = 0;%rand(1);
ds1 = ds;
ds2 = ds;
s1 = r + hcm + ds1;
s2 = r + hcm + ds2;
% construct rotation matrix for angle t relative to [0;0;1] in the yz plane
vx = [1;0;0];
vy = [0;cos(t);sin(t)];
vz = [0;-sin(t);cos(t)];
Rt = [vx vy vz];
% position of the contact points
po1 = [0;0;s1-r]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;s2-r]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;
% contact normals
n1 = [-sin(ph1),-cos(ph1)*cos(t),-cos(ph1)*sin(t)].';
n2 = [-sin(ph2),-cos(ph2)*cos(t),-cos(ph2)*sin(t)].';
n3 = [0,-sin(t),cos(t)].';

com = Rt*com;


%% force magnitudes
N = [n1 n2 n3];
g = [0 0 -1].';
L = -N\g;

%% finding p3 for coplanar p1p2
%  syms ds h3
if ds1 == ds2 &&ph1==-ph2
    h12 = cos(t)*ds;
    h3 = ds*(cos(t)^2+1)/cos(t);
    s3 = (h3-h12)*sin(t);
    if s3>=0 && s3<=r
        p3 = s3*[0;cos(t);sin(t)];
    else
        p3 = [];
        continue
    end
else
    A = [n1(1:2) -n2(1:2)];
    dp = p2-p1;
    sol = A\dp(1:2);
    inter1 = p1+sol(1)*n1;
    % inter2 = p2+sol(2)*n2;
    %     h1 = inter1(3)-com(3);
    %     h2 = inter2(3)-com(3);
    lxy = inter1(1:2)-com(1:2);
    h1 = -(ds1*cos(ph1)*sin(ph2) - ds2*cos(ph1)*sin(ph2) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
    h2 = -(ds1*cos(ph2)*sin(ph1) - ds2*cos(ph2)*sin(ph1) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
    h3 = -(ds1*cos(ph1)*sin(ph2) - ds2*cos(ph2)*sin(ph1) - ds1*cos(ph2)*cos(t)^2*sin(ph1) + ds2*cos(ph1)*cos(t)^2*sin(ph2))/(sin(ph1 - ph2)*cos(t));
    inter3 = [lxy;h3]+com;
    dist3 = n3.'*inter3;
    p3 = inter3-n3*dist3;
end
 
%% plotting the contacts and force lines
plotview = 1;
if plotview
figure();
subplot(1,2,1) % yz plane
hold on

%plot rectangle view of the cylinder
v1 = -r*[cos(t);sin(t)];
v2 = r*[cos(t);sin(t)];
v3 = r*[cos(t);sin(t)]+h*[-sin(t);cos(t)];
v4 = -r*[cos(t);sin(t)]+h*[-sin(t);cos(t)];
vertex = [v1 v2 v3 v4 v1];
plot(vertex(1,:),vertex(2,:),'b');
% plot axis of cylinder
axisends = [v1+v2 (v3+v4)/2];
plot(axisends(1,:),axisends(2,:),'b--');
% mark com
plot(com(2),com(3),'*r');
plot([com(2),com(2)],[0,com(3)*2],'r--');
%plot contact points
plot(p1(2),p1(3),'ok','markerSize',4,'lineWidth',4);
plot(p2(2),p2(3),'ok','markerSize',4,'lineWidth',4);
plot(p3(2),p3(3),'ok','markerSize',4,'lineWidth',4);
% plot force lines
f1 = [p1 p1+n1];
plot(f1(2,:),f1(3,:),'k--');
f2 = [p2 p2+n2];
plot(f2(2,:),f2(3,:),'k--');
f3 = [p3 p3+2*n3];
plot(f3(2,:),f3(3,:),'k--');
xlabel('Y')
ylabel('Z')

axis equal
hold off

%
subplot(1,2,2) % yx plane
hold on

% plot axis of cylinder
e1 = [0;0];
e2 = [-h*sin(t);0];
axisends = [e1 e2];
plot(axisends(1,:),axisends(2,:),'b--');

% mark com
plot(com(2),com(1),'*r')

%plot contact points
plot(p1(2),p1(1),'ok','markerSize',4,'lineWidth',4)
plot(p2(2),p2(1),'ok','markerSize',4,'lineWidth',4)
plot(p3(2),p3(1),'ok','markerSize',4,'lineWidth',4)

% plot force lines
f1 = [p1 p1+n1];
plot(f1(2,:),f1(1,:),'k--');
f2 = [p2 p2+n2];
plot(f2(2,:),f2(1,:),'k--');
f3 = [p3 p3+2*n3];
plot(f3(2,:),f3(1,:),'k--');
xlabel('Y')
ylabel('X')
hold off
axis equal
end
%% checking stability
[flag3, eigv] = MasonConditionFunction(r,hcm,t,ph,ds)
[flag1,val] = MinCondition2D(hcm,r,ph,t,ds)
if ~flag1 && flag3  %1 % flag3 %sum(eigv>0)>1 ||~flag1
    draw2D(r,t,h,com,p1,p2,p3,n1,n2,n3)
    disp('need to draw 3-contact plane, or compute wrenches for object at slightly changed orientation')
    %% draw 3-contact plane
    base = norm(p1-p2);
    base_middle = (p1+p2)/2;
    v2p3 = base_middle-p3;
    height = norm(v2p3);
    slope = v2p3(3)/(v2p3(1)^2+v2p3(2)^2)^0.5;
    dph = 0.01; dth = 0.01;
    res = 4;
    ph1_vector = ph1-res*dph:dph:ph1+res*dph;
    th1_vector = t-res*dth:dth:t+res*dth;
    eqgrade = zeros(length(ph1_vector),length(th1_vector));
    ComDepth = zeros(length(ph1_vector),length(th1_vector));
    RotAngle = zeros(length(ph1_vector),length(th1_vector));
    for k=1:length(ph1_vector)
        ph1t = ph1_vector(k);
        for j=1:length(th1_vector)
            tt = th1_vector(j);
            [p1c,ph1c,p2c,ph2c,p3c,hc] = TripleContactStance(ph1t,tt,base,height,slope,r,hcm);
            if p3c~=inf
                ComDepth(k,j) = p3(3)+hc;
                p3e = equilibriumStance(p1c,p2c,ph1c,ph2c,tt,hcm);
                eqgrade(k,j) = norm(p3e-p3c);
                p2p1 = p2c-p1c;
                RotAngle(k,j) = wrapToPi(-atan2(p2p1(2),p2p1(1)));
            else
                ComDepth(k,j) = inf;
                eqgrade(k,j) = inf;
            end
            
        end
    end
    %%
    f = figure;
    contour(ph1_vector,th1_vector,ComDepth.',10); %[linspace(min(min(ComDepth)),max(max(ComDepth(ComDepth<inf))),50),min(min(ComDepth))+logspace(-4,-1,50)])
    % contour(ph1_vector,th1_vector,ComDepth.',100);
    xlabel('\phi_1')
    ylabel('\theta_1')
    zlabel('COM Height')
    cb=colorbar;
    set(get(cb,'label'),'string','U','interpreter','latex','FontSize',14,'rotation',0);
    x = 5
end

end

