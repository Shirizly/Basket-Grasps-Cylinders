function drawCliff3D(cliff_s3,cliff_phi3,cliff_height,col)
Ni = length(cliff_s3);
cliff_s3_p = [cliff_s3 cliff_s3(Ni:-1:1)];
cliff_phi3_p = [cliff_phi3 -cliff_phi3(Ni:-1:1)];

cliff_x = cliff_s3_p.*sin(cliff_phi3_p);
cliff_y = cliff_s3_p.*cos(cliff_phi3_p);
hold on
if ~isempty(cliff_height)
    cliff_height = [cliff_height cliff_height(Ni:-1:1)];
    plot3(cliff_x,cliff_y,cliff_height,col,'lineWidth',4)
else
    plot(cliff_x,cliff_y,col,'lineWidth',4)
end
axis equal
end