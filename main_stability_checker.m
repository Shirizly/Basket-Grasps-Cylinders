num = 0;
curv = 1;
suf = 1;
N = 10; Nt1 = 50; Np1 = 50; Nds = Nt1;
r = 0.5;
hcm = 1;
% t1_range = linspace(1.151,1.323,Nt1);
t1_range = linspace(0.01,1.55, Nt1);%
% t1_range = linspace(pi()/20,9*pi()/20,Nt1);
% t1 = 0.48;
phi_range = linspace(pi()*(0.5+0.5*pi()/Np1),pi()*(1-0.5*pi()/Np1),Np1);
% ds_range = linspace(0.1*hcm,0.3*hcm,Nds);
escape_stance_phi3 = -inf*ones(length(t1_range),length(phi_range));
flagNumA = zeros(Nt1,Np1);
flagCurvA = zeros(Nt1,Np1);
flagSuffA = zeros(Nt1,Np1);
valNumA = zeros(Nt1,Np1);
valCurvA = zeros(Nt1,Np1);
val2DA = zeros(Nt1,Np1);
escape_height_arrayDSRec = zeros(Nt1,Np1);
escape_height_arrayDSCirc = zeros(Nt1,Np1);
escape_height_arraySS = zeros(Nt1,Np1);
i=1;
TSNonTrivialCounter = 0;
dsv = -0.51;%linspace(-hcm/2,-0.01,2);
% dsv = [-0.2,-0.1,-0.05,-0.03,-0.01];
maxi = 0;
for r=1%linspace(0.1,50,10)
for ds = dsv
tic
ds
flagNumA = -inf*ones(Nt1,Np1);
flagCurvA = -inf*ones(Nt1,Np1);
flagSuffA = -inf*ones(Nt1,Np1);
valNumA = -inf*ones(Nt1,Np1);
valCurvA = -inf*ones(Nt1,Np1);
% val2DA = zeros(Nt1,Np1);
for j = 1:length(t1_range)
    t1 = t1_range(j); %#ok<NOPTS>
%     f(i,j) = figure; %#ok<SAGROW>
%     hold on
%     ds = ds_range(j) %#ok<NOPTS>
    %     ds = 0.1*hcm;
    com_0 = [0;0;hcm];
    
    for k=1:length(phi_range)
        phi = phi_range(end+1-k);
        %             [flag3, eigv] = MasonConditionFunction(r,hcm,tt,phi,ds);
        %             if ~flag3
        %                 continue
        %             end
        ph1 = phi; ph2 = -phi;
        vx = [1;0;0];
        vy = [0;cos(t1);sin(t1)];
        vz = [0;-sin(t1);cos(t1)];
        Rt = [vx vy vz];
        
        com = Rt*com_0;
%         com_check = GraspComPos(ds,t1,r,phi);
%         norm(com-com_check)
        % position of the contact points
        po1 = [0;0;hcm+ds]+r*[sin(ph1);cos(ph1);0];
        p1 = Rt*po1;
        po2 = [0;0;hcm+ds]+r*[sin(ph2);cos(ph2);0];
        p2 = Rt*po2;
%         p3 = equilibriumStance(p1,p2,ph1,ph2,t1,hcm);
        p3 = ds*tan(t1)*[0;cos(t1);sin(t1)];
%         check = 0.5*(p1+p2-2*(hcm+ds)*com)-p3;
%         check(2)
        if norm(p3)<r
%             [flag2D,val2D] = MinCondition2D(hcm,r,phi,t1,ds,p1,p2,p3);
            [Curvature,val_curv] = curvatureConFast(t1,phi,ds,r,hcm);
            Numerical = 0;
%             [Numerical,Numval] = NumericalStabilityCheck(r,p1,p2,p3,t1,N);
            boundval1 = -hcm/(1+sqrt(1+r^2));%(1+sqrt(1+ds*tan(t1)));
            boundval2 = r*cos(phi)/(2*tan(t1));
            if ds>max(boundval1,boundval2)
                flagSuff=1;
            else
                flagSuff=0;
            end
%             if (Numerical&&Curvature)~=(Numerical||Curvature)
%                   disp('disagreement')
%             end
            if flagSuff==1 && Curvature==0
                disp('Error: sufficient condition insufficient');
            end
            boundval3 = -hcm/2;
            if ds<min(boundval2,boundval3) && Curvature==1
                disp('not sufficient for instability');
            end
            
            flagNumA(j,end+1-k) = Numerical;
            flagCurvA(j,end+1-k) = Curvature;
            flagSuffA(j,end+1-k) = flagSuff;
%             valNumA(j,end+1-k) = Numval;
%             valCurvA(j,end+1-k) = max(val_curv);
%             val2DA(j,k) = -val2D;            
            
        end
    end
    
end
toc
%%




dsc = 0;
for i=1:dsc
subFolderName = '\Depth_r_05_h_2_ds\zoom';
fig1 = figure;
[~,cf] = contourf(ds_range,phi_range,escape_height_arrayDSRec.',30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\delta s$','interpreter','latex','FontSize',20)
ylabel('$\phi_1$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U_{2c,1}$','interpreter','latex','FontSize',20,'rotation',0);
xlim([min(ds_range),max(ds_range(escape_height_arraySS(:,1)>0))])
savefig(fig1,[pwd '\' subFolderName '\2c,1'])
saveas(fig1,[pwd '\' subFolderName '\2c,1.eps'],'epsc')
saveas(fig1,[pwd '\' subFolderName '\2c,1.jpeg'],'jpeg')

fig2 = figure
[~,cf] = contourf(ds_range,phi_range,escape_height_arrayDSCirc.',30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\delta s$','interpreter','latex','FontSize',20)
ylabel('$\phi_1$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U_{2c,2}$','interpreter','latex','FontSize',20,'rotation',0);
xlim([min(ds_range),max(ds_range(escape_height_arraySS(:,1)>0))])
savefig(fig2,[pwd '\' subFolderName '\2c,2'])
saveas(fig2,[pwd '\' subFolderName '\2c,2.eps'],'epsc')
saveas(fig2,[pwd '\' subFolderName '\2c,2.jpeg'],'jpeg')

fig3 = figure;
[~,cf] = contourf(ds_range,phi_range,escape_height_arraySS.',50);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\delta s$','interpreter','latex','FontSize',20)
ylabel('$\phi_1$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U_{1c}$','interpreter','latex','FontSize',20,'rotation',0);
xlim([min(ds_range),max(ds_range(escape_height_arraySS(:,1)>0))])
savefig(fig3,[pwd '\' subFolderName '\1c'])
saveas(fig3,[pwd '\' subFolderName '\1c.eps'],'epsc')
saveas(fig3,[pwd '\' subFolderName '\1c.jpeg'],'jpeg')

fig4 = figure;
[~,cf] = contourf(ds_range,phi_range,escape_height_arrayTS.',30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\delta s$','interpreter','latex','FontSize',20)
ylabel('$\phi_1$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U_{3c}$','interpreter','latex','FontSize',20,'rotation',0);
xlim([min(ds_range),max(ds_range(escape_height_arraySS(:,1)>0))])
savefig(fig4,[pwd '\' subFolderName '\3c'])
saveas(fig4,[pwd '\' subFolderName '\3c.eps'],'epsc')
saveas(fig4,[pwd '\' subFolderName '\3c.jpeg'],'jpeg')

fig5 = figure;
[~,cf] = contourf(ds_range,phi_range,escape_height_array.',30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\delta s$','interpreter','latex','FontSize',20)
ylabel('$\phi_1$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U$','interpreter','latex','FontSize',20,'rotation',0);
xlim([min(ds_range),max(ds_range(escape_height_arraySS(:,1)>0))])
savefig(fig5,[pwd '\' subFolderName '\U'])
saveas(fig5,[pwd '\' subFolderName '\U.eps'],'epsc')
saveas(fig5,[pwd '\' subFolderName '\U.jpeg'],'jpeg')
end
%% plot and save figures of phi-theta
theta = 1;
for i = 1:theta
    
%     subFolderName = [pwd '\Depth_r_0' num2str(r*10) '_ds_0' num2str(ds*100)];
%     mkdir(char(subFolderName))
if num
fig1 = figure;
% hold on
% [~,~] = contourf(t1_range,phi_range,valNumA.',[0]);
[~,cf] = contourf(t1_range,phi_range,flagNumA.',30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\theta$','interpreter','latex','FontSize',20)
ylabel('$\phi$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','Num','interpreter','latex','FontSize',20,'rotation',0);
% xlim([min(t1_range),max(t1_range(valNumA(:,1)~=0))])
title(['distance from COM along the axis = ' num2str(ds)])
% savefig(fig1,[pwd '\' subFolderName '\2c,1'])
% saveas(fig1,[pwd '\' subFolderName '\2c,1.eps'],'epsc')
% saveas(fig1,[pwd '\' subFolderName '\2c,1.png'],'png')
end
if curv
fig2 = figure;
% hold on
% [~,~] = contourf(t1_range,phi_range,valCurvA.',[0]);
[~,cf] = contourf(t1_range,phi_range,flagCurvA.',30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\theta$','interpreter','latex','FontSize',20)
ylabel('$\phi$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','Curv','interpreter','latex','FontSize',20,'rotation',0);
% xlim([min(t1_range),max(t1_range(valCurvA(:,1)~=0))])
title(['distance from COM along the axis = ' num2str(ds)])
% savefig(fig2,[pwd '\' subFolderName '\2c,2'])
% saveas(fig2,[pwd '\' subFolderName '\2c,2.eps'],'epsc')
% saveas(fig2,[pwd '\' subFolderName '\2c,2.png'],'png')
end
if suf
fig3 = figure;
[~,cf] = contourf(t1_range,phi_range,flag2DA.',50);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\theta$','interpreter','latex','FontSize',20)
ylabel('$\phi$','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','suff','interpreter','latex','FontSize',20,'rotation',0);
% xlim([min(t1_range),max(t1_range(val2DA(:,1)~=0))])
title(['distance from COM along the axis = ' num2str(ds)])
% savefig(fig3,[pwd '\' subFolderName '\1c'])
% saveas(fig3,[pwd '\' subFolderName '\1c.eps'],'epsc')
% saveas(fig3,[pwd '\' subFolderName '\1c.png'],'png')
end
% fig4 = figure;
% [~,cf] = contourf(t1_range,phi_range,escape_height_arrayTS.',30);
% % set(cf,'LineColor','none')
% cb = colorbar;
% xlabel('$\theta_1$','interpreter','latex','FontSize',20)
% ylabel('$\phi_1$','interpreter','latex','FontSize',20)
% set(get(cb,'label'),'string','$\:U_{3c}$','interpreter','latex','FontSize',20,'rotation',0);
% xlim([min(t1_range),max(t1_range(escape_height_arraySS(:,1)>0))])
% title(['distance from COM along the axis = ' num2str(ds)])
% % savefig(fig4,[pwd '\' subFolderName '\3c'])
% % saveas(fig4,[pwd '\' subFolderName '\3c.eps'],'epsc')
% % saveas(fig4,[pwd '\' subFolderName '\3c.png'],'png')
% 
% fig5 = figure;
% [~,cf] = contourf(t1_range,phi_range,escape_height_array.',30);
% % set(cf,'LineColor','none')
% cb = colorbar;
% xlabel('$\theta_1$','interpreter','latex','FontSize',20)
% ylabel('$\phi_1$','interpreter','latex','FontSize',20)
% set(get(cb,'label'),'string','$\:U$','interpreter','latex','FontSize',20,'rotation',0);
% xlim([min(t1_range),max(t1_range(escape_height_arraySS(:,1)>0))])
% title(['distance from COM along the axis = ' num2str(ds)])
% % savefig(fig5,[pwd '\' subFolderName '\U'])
% % saveas(fig5,[pwd '\' subFolderName '\U.eps'],'epsc')
% % saveas(fig5,[pwd '\' subFolderName '\U.png'],'png')
end
end
end

function [phi3_range,height,centersolA,comsolA] = cutSectionOffCliffCurve(phi3_range,height,centersolA,comsolA,min_ind,minNonTSHeight)
if isempty(minNonTSHeight)
    numer_deriv = diff(height);
    if min_ind ==1
        cutoff1 = 1;
    else
    cutoff1 = find(numer_deriv(1:min_ind-1)>0,1,'last');
    end
    if isempty(cutoff1)
        cutoff1 = 1;
    end
    if min_ind==length(height)
        cutoff2 = length(height);
    else
    cutoff2 = find(numer_deriv(min_ind:end)<0,1,'first');
    end
    if isempty(cutoff2)
        cutoff2 = length(height);
    end
end
phi3_range(cutoff1:cutoff2) = [];
height(cutoff1:cutoff2) = [];
centersolA(cutoff1:cutoff2,:) = [];
comsolA(cutoff1:cutoff2,:) = [];
end