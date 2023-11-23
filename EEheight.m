function sol = EEheight(dp,nj,nk,vj,vk,dels1,dels2)
% Function returns the obect configuration at Edge-Edge maximum height equilibrium
% stances. 
dp1 = dp(1);
dp2 = dp(2);
nj1 = nj(1);
nj2 = nj(2);
nk1 = nk(1);
nk2 = nk(2);
vj1 = vj(1);
vj2 = vj(2);
vk1 = vk(1);
vk2 = vk(2);
% transform the quadratic equation in cosine and sine to 4th order poly.
Avector = EEheightAvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
a1 = Avector(1);
a2 = Avector(2);
a3 = Avector(3);
a4 = Avector(4);
a5 = Avector(5);
% find the coefficients for the canonical form of 4th order poly.
p1 = a1 - a3;
p2 = 2*a5 - 2*a2;
p3 = 4*a4 - 2*a1;
p4 = 2*a2 + 2*a5;
p5 = a1 + a3;
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
eps = 1E-7;
theta = ones(4,1);
for i=1:4
    theta(i) = NaN;
    if abs(imag(x(i)))<eps
        theta(i) = 2*atan2(real(x(i)),1);
    end
end

% calculate the second derivative of the height according to theta
Cvector = EEheightCvector(dp1,dp2,nj1,nj2,nk1,nk2,vj1,vj2,vk1,vk2);
Terms = zeros(4,5);
diff2d = zeros(4,1);
diffd = zeros(4,1);
for i=1:4
Terms(i,:) = [ cos(theta(i))^2, cos(theta(i))*sin(theta(i)), cos(theta(i)), sin(theta(i))^2, sin(theta(i))];

diff2d(i) = Terms(i,:)*Cvector.';

diffd(i) = Terms(i,:)*Avector.';
end
if any(abs(diffd)>1E-5)
    disp('error in calculating equilibrium');
end
%% find if the contact is within the actual edge and check feasibility
% also find the COM position relative to finger 1

e = [0,1];
J = [0 1; -1 0];
A = [-J*nj J*nk];
B = [J*nj J*nk];
dv = vj-vk;
sv = vj+vk;
s = zeros(2,4);
d = zeros(2,4);
feasicond = zeros(4,1);
for i=1:4
    R = [cos(theta(i)) -sin(theta(i));sin(theta(i)) cos(theta(i))];
    feasicond(i) = ((e*R*nj>0)||(e*R*nk>0))&&(([1,0]*R*nj)*([1,0]*J*R*nk)<0);
    if ~any(isnan(R))
        lamda = [R*nj R*nk]\[0;1];
        feasicond(i) = feasicond(i)&& all(lamda>0);
    end
    s(:,i) = A\(R.'*dp+dv);
    d(:,i) = 0.5*(dp-R*(sv+B*s(:,i)));
    % remove solutions that don't fit the two criterions:
    if  any(s(:,i)<0) || s(1,i)>dels1 || s(2,i)>dels2 || diff2d(i)<0 || ~feasicond(i)
        theta(i) = NaN;
    end
end

sol = [];
% prepare the solution vectors for output
for i=1:4
    if ~isnan(theta(i))
        sol = [sol,[d(:,i);theta(i)]];
    end
end

end
