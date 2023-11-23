function drawEllipse(a,b,c,theta,p1,p2)
hold on
plot([p1(1) p2(1)],[p1(2) p2(2)],'.r')
t = 0:0.001:2*pi()+0.001;
R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
el = R*[a*cos(t);b*sin(t)]+c;
plot(el(1,:),el(2,:),'b');
plot(c(1),c(2),'+b');
axis equal
grid on
end