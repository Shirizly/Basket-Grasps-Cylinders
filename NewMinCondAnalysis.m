a = 1.5;
b = 1;

p1 = [1.5;0.708];
p2 = [1.5;-0.708];
figure 
th1 = 0;
% tic
c1 = findEllipseCenter(a,b,th1,p1,p2);
% toc
drawEllipse(a,b,c1,th1,p1,p2)
th2 = 0.1;
% tic
c2 = findEllipseCenter(a,b,th2,p1,p2);
drawEllipse(a,b,c2,-th2,p1,p2)
% toc
pvec = p2-c1;
pvecCircle = [1/a 0;0 1/b]*pvec;
t = atan2(pvecCircle(2),pvecCircle(1));
CondAlon2 = c2(1)-c1(1)
K = a*b/((a^2*sin(t)^2+b^2*cos(t)^2)^3/2)

%% run over slice of configuration space
r = 0.5;
h = 2;
hcm = 0.5*h;
theta_range = 0.01:0.01:0.45;%pi()/2-0.1;
phi_range = pi()/2+0.02:0.02:pi()-0.02;
figure
hold on
for th = theta_range
    CondAlon2 = [];
    for ph = phi_range        
        ph1 = ph; ph2 = -ph;
        com = [0;0;hcm];
        ds = 0.1*h;
        ds1 = ds;
        ds2 = ds;
        s1 = r + hcm + ds1;
        s2 = r + hcm + ds2;
        % construct rotation matrix for angle t relative to [0;0;1] in the yz plane
        vx = [1;0;0];
        vy = [0;cos(th);sin(th)];
        vz = [0;-sin(th);cos(th)];
        Rt = [vx vy vz];
        % position of the contact points
        po1 = [0;0;s1-r]+r*[sin(ph1);cos(ph1);0];
        p1 = Rt*po1;
        po2 = [0;0;s2-r]+r*[sin(ph2);cos(ph2);0];
        p2 = Rt*po2;
        h12 = cos(th)*ds;
        h3 = ds*(cos(th)^2+1)/cos(th);
        s3 = (h3-h12)*sin(th);
        if s3>0 && s3<r
            p3 = s3*[0;cos(th);sin(th)];
        else
            p3 = [];
            CondAlon2(end+1) = inf;
            continue
        end
        % ellipse analysis
        p1e = p1(1:2); p2e = p2(1:2); 
        a = r/cos(th);
        b = r;
        c1 = findEllipseCenter(a,b,0,p1e,p2e);
        c2 = findEllipseCenter(a,b,0.1,p1e,p2e);
        CondAlon2(end+1) = (c2(2)-c1(2));
        
        
    end
    plot(phi_range,CondAlon2,'DisplayName',['theta=' num2str(round(th,2))])
end
xlabel('\phi_1,2')
ylabel('stability condition')
title('stability condition for various values of \theta_1')
