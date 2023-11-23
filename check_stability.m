function flag = check_stability(r,p1,p2,p3,t,hcm)
flag = 0;
s3 = norm(p3);
h = CylinderComPos(r,hcm,p1,p2,p3,t,0,s3);
N=10;
delta = s3/10;
minht = inf;
for i=1:N
    s3t = s3 + delta*cos(2*pi()*i/N);
    phi3t = atan2(delta*sin(2*pi()*i/N),s3);
    ht = CylinderComPos(r,hcm,p1,p2,p3,t,phi3t,s3t);
    minht = min(minht,ht);
     h-minht
end
if h<minht
    flag = 1;
end
end