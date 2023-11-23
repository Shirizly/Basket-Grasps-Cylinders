phi = 3.0916;
hcm = 0.8;
tt = pi()/10;
ds = 0.4*hcm;
%             [flag3, eigv] = MasonConditionFunction(r,hcm,tt,phi,ds);
%             if ~flag3
%                 continue
%             end
ph1 = phi; ph2 = -phi;
vx = [1;0;0];
vy = [0;cos(tt);sin(tt)];
vz = [0;-sin(tt);cos(tt)];
Rt = [vx vy vz];

com = Rt*com;
% position of the contact points
po1 = [0;0;hcm+ds]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;hcm+ds]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;
p3 = equilibriumStance(p1,p2,ph1,ph2,tt,hcm);
if norm(p3)<r
    [phi3_range,height,theta1,theta2] = CylinderCliff(r,hcm,p1,p2,p3,tt,N);
    height = height - com(3);
    [min_height,min_ind] = min(height);
    if abs(phi3_range(min_ind))>1E-1
        %%
        [minCOM,mintheta1,mintheta2] = CylinderCliffPoint(r,hcm,p1,p2,p3,tt,N,0,[theta1(101),theta2(101)]);
%         Draw3D(r,hcm*2,hcm,tt,p1,p2,p3,ph1,ph2)
        Draw3DEscape(r,hcm*2,hcm,minCOM,mintheta1,mintheta2,p1,p2,p3)
        disp('non-trivial escape')
        %%
    end
    if min_height>0
        escape_stance_phi3(j,k) = phi3_range(height==min(height));
        plot(phi3_range,height,'.','DisplayName',['theta1 ' num2str(tt) ' phi ' num2str(phi)])
    end
end