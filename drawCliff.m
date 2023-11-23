function drawCliff(cliff_s3,cliff_phi3,col)
Ni = length(cliff_s3);
cliff_s3_p = [cliff_s3 cliff_s3(Ni:-1:1)];
cliff_phi3_p = [cliff_phi3 -cliff_phi3(Ni:-1:1)];

cliff_x = cliff_s3_p.*sin(cliff_phi3_p);
cliff_y = cliff_s3_p.*cos(cliff_phi3_p);
hold on
plot(cliff_x,cliff_y,'color',col,'lineWidth',4)


hold off
axis equal
end