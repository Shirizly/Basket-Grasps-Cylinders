function fig = drawCliff3DApprox(r,p1,p2,p3,centerA,linestyle,fig)
if isempty(fig)
    fig = draw3DSpheres(r,p1,p2,p3);
end
hold on
plot3(centerA(:,1),centerA(:,2),centerA(:,3),linestyle,'lineWidth',4,'markerSize',2);
% for i=1:size(nsolA,1)
%     centerA = centersolA(i,:).';
%     comA = comsolA(i,:).';
%     axisA = [centerA comA];
%     
%     plot3(axisA(1,:),axisA(2,:),axisA(3,:),'.k','lineWidth',4,'markerSize',2);
% end
hold off
end