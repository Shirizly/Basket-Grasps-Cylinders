function sol = EEheight_new(r,hcm,p_c1,p_c2,dels1,dels2)
% Function returns the obect configuration at Edge-Edge maximum height equilibrium
% stances. 
debug = 0;
p_c11 = p_c1(1);
p_c12 = p_c1(2);
p_c21 = p_c2(1);
p_c22 = p_c2(2);

% find the coefficients for the canonical form of 4th order poly.
poly4 = [p_c21*1i - p_c11*1i - p_c12 + p_c22 , - r*1i + hcm , 0 , - r*1i - hcm , p_c21*1i - p_c22 - p_c11*1i + p_c12];
p1 = poly4(1);
p2 = poly4(2);
p3 = poly4(3);
p4 = poly4(4);
p5 = poly4(5);
% solve the 4th order poly.
p = (8*p1*p3-3*p2^2)/(8*p1^2);

q = (p2^3-4*p1*p2*p3+8*p1^2*p4)/(8*p1^3);

del0 = p3^2-3*p2*p4+12*p1*p5;

del1 = 2*p3^3-9*p2*p3*p4+27*p2^2*p5+27*p1*p4^2-72*p1*p3*p5;

Q = ((del1+sqrt(del1^2-4*del0^3))/2)^(1/3);

S = 0.5*sqrt(-2*p/3+1/(3*p1)*(Q+del0/Q));


x(1) = -p2/(4*p1)-S+0.5*sqrt(-4*S^2-2*p+q/S);
x(2) = -p2/(4*p1)-S-0.5*sqrt(-4*S^2-2*p+q/S);
x(3) = -p2/(4*p1)+S+0.5*sqrt(-4*S^2-2*p-q/S);
x(4) = -p2/(4*p1)+S-0.5*sqrt(-4*S^2-2*p-q/S);
% check each solution for realness, and if real, calculate theta
eps = 1E-9;
theta = ones(4,1);
for i=1:4
        theta(i) = - log(x(i))*1i;
        if imag(theta(i))>eps
            theta(i) = inf;
        else
            theta(i) = wrapToPi(real(theta(i)));
        end
            
end

% calculate the second derivative of the height according to theta
dpdtheta = zeros(4,1);
d2pdtheta2 = zeros(4,1);
% for i=1:4
% dpdtheta(i) = p_c11*cos(2*theta(i)) - p_c21*cos(2*theta(i)) + p_c12*sin(2*theta(i)) - p_c22*sin(2*theta(i)) + r*cos(theta(i)) - hcm*sin(theta(i));
% d2pdtheta2(i) = 2*p_c12*cos(2*theta(i)) - 2*p_c22*cos(2*theta(i)) - 2*p_c11*sin(2*theta(i)) + 2*p_c21*sin(2*theta(i)) - hcm*cos(theta(i)) - r*sin(theta(i));
% end
% if any(abs(dpdtheta)>1E-7)
%     disp('error in calculating equilibrium');
% end

    
%% find if the contact is within the actual edge and check feasibility
% also find the COM position relative to finger 1

e = [0,1];
J = [0 1; -1 0];
% A = [-J*nj J*nk];
% B = [J*nj J*nk];
% dv = vj-vk;
% sv = vj+vk;
nj = [0;1];
nk = [1;0];
s = zeros(2,4);
d = zeros(2,4);
feasicond = zeros(4,1);
for i=1:4
    tsol = theta(i);
    dpdtheta(i) = p_c11*cos(2*tsol) - p_c21*cos(2*tsol) + p_c12*sin(2*tsol) - p_c22*sin(2*tsol) + r*cos(tsol) - hcm*sin(tsol);
    d2pdtheta2(i) = 2*p_c12*cos(2*tsol) - 2*p_c22*cos(2*tsol) - 2*p_c11*sin(2*tsol) + 2*p_c21*sin(2*tsol) - hcm*cos(tsol) - r*sin(tsol);
    if (abs(dpdtheta(i))>1E-7)
        disp('error in calculating equilibrium');
    end
    R = [cos(theta(i)) -sin(theta(i));sin(theta(i)) cos(theta(i))];
    feasicond(i) = ((e*R*nj>0)||(e*R*nk>0))&&(([1,0]*R*nj)*([1,0]*J*R*nk)<0); %reaction forces span the origin
    if ~any(isnan(R))
        lambda = [R*nj R*nk]\[0;1]; %magnitudes of reaction forces
        feasicond(i) = feasicond(i)&& all(lambda>0); %make sure reaction forces are pushing only
    end
    s(:,i) = [2*r + p_c11*cos(tsol) - p_c21*cos(tsol) + p_c12*sin(tsol) - p_c22*sin(tsol);...
        p_c12*cos(tsol) - p_c22*cos(tsol) - p_c11*sin(tsol) + p_c21*sin(tsol)];
    d(:,i) = [p_c21 - hcm*sin(tsol) + cos(tsol)*(r + p_c11*cos(tsol) - p_c21*cos(tsol) + p_c12*sin(tsol) - p_c22*sin(tsol));...
        p_c22 + hcm*cos(tsol) + sin(tsol)*(r + p_c11*cos(tsol) - p_c21*cos(tsol) + p_c12*sin(tsol) - p_c22*sin(tsol))];
    if debug
        figure;hold on;
        plot(p_c11,p_c12,'bo')
        plot(p_c21,p_c22,'bo')
        edge1 = [p_c2+s(1,i)*R*J*nj,p_c2-(dels1-s(1,i))*R*J*nj];
        edge2 = [p_c1+s(2,i)*R*J*nk,p_c1-(dels2-s(2,i))*R*J*nk];
        plot(edge1(1,:),edge1(2,:),'k');
        plot(edge2(1,:),edge2(2,:),'k');
        plot(d(1,i),d(2,i),'+r');
        axis equal
    end
    % remove solutions that don't fit the two criterions:
    if  any(s(:,i)<0) || s(1,i)>dels1 || s(2,i)>dels2 || d2pdtheta2(i)<0 || ~feasicond(i)
        theta(i) = NaN;
    else
        break
    end
end

sol = [];
% prepare the solution vectors for output
for i=1:4
    if ~isnan(theta(i))
        sol = [sol,[d(:,i);theta(i)]];
        break
    end
end

end
