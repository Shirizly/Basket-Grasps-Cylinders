tic
N = 50; Nt1 = 10; Np1 = 10; Nds = Nt1;
rv = 5;
r = rv;
hcm = 1;
t1_rangev = linspace(0,1.55, Nt1);%
% t1_range = linspace(pi()/20,9*pi()/20,Nt1);
% t1 = 0.48;
phi_range = linspace(pi()*21/40,pi()*39/40,Np1);%
phi = 2;
% ds_range = linspace(0.1*hcm,0.3*hcm,Nds);
dsv = linspace(0,0.5,Nds);
ds = 0;
maxi = 0;
count_no_circ = 0;
tolerance = 1E-4;
% for 

    t1_range = t1_rangev;
%     t1_range(t1_rangev>atan(r/ds))=[];
    Nt1a = length(t1_range);
    Np1a = length(phi_range);
    escape_stance_phi3 = -inf*ones(Nt1a,Np1a);
    escape_height_array = zeros(Nt1a,Np1a);
    escape_height_arrayTSA = zeros(Nt1a,Np1);
    escape_height_arrayTS = zeros(Nt1a,Np1a);
    escape_height_arrayDSRec = zeros(Nt1a,Np1a);
    escape_height_arrayDSCirc = zeros(Nt1a,Np1a);
    escape_height_arraySS = zeros(Nt1a,Np1a);
    escape_type = zeros(Nt1a,Np1a);
    i=1;
    TSNonTrivialCounter = 0;
    count_numer = 0;
    count_approx = 0;
    for j = 1:length(t1_range)
        %%
        t1 = t1_range(j) %#ok<NOPTS>
        %     f(i,j) = figure; %#ok<SAGROW>
        %     hold on
        %     ds = ds_range(j) %#ok<NOPTS>
        %     ds = 0.1*hcm;
        
        
%         for k=1:length(dsv)
%             ds = dsv(k);
%             tic
        for k=1:length(phi_range)
            phi = phi_range(k);
            
            p3 = ds*tan(t1)*[0;cos(t1);sin(t1)];
            if norm(p3)<r %necessary for feasibility of grasp, otherwise base finger needs to support object at a point that doesn't exist
                [p1,p2,com] = grasp_details(r,hcm,phi,t1,ds);
                P = [p1,p2,p3];
                [flag_curv,val_curv] = curvatureConFast(t1,phi,ds,r,hcm); %analyitcally derived stability condition
                %              [flag,val] = NumericalStabilityCheck(r,p1,p2,p3,t1,N); %numerical condition, slow
                if flag_curv
                    %% single-support escape stances
                    min_heightSS = singleSupportEscape(r,hcm,p1,p2,p3);
                    min_depthSS = min_heightSS - com(3);
                    
                    %% double-support escape stances
%                     Draw3DGrasp(r,hcm*2,hcm,t1,p1,p2,p3);
                    min_heightDS = doubleSupportEscape(r,hcm,p1,p2,p3);
                    min_depthDS = min_heightDS - com(3);
                    if p3(3)>p1(3) && abs(p3(2)-p1(2))<hcm
                        count_no_circ = count_no_circ+1;
                    end
                    No_skip = 1;
                    if No_skip
                    %% triple-support escape stances
                    if norm(p1-p3)>2*r
                    [phi3_range,height,nsolA,comsolA,centersolA] = CylinderCliff3(r,hcm,p1,p2,p3,t1,phi,N);
                    for i=1:size(centersolA,1)
                        s3_range(i) = norm(centersolA(i,:));
                    end
                    count_approx = count_approx+1;
                    else
                        [phi3_range,s3_range,height,comsolA,nsolA] = CylinderCliffNumericalBinary(r,hcm,t1,P,N,tolerance,0);
                        %                     [phi3_range,s3_range,height] = CylinderTS2(r,hcm,P,t1,phi,50,1);
                        count_numer = count_numer+1;
                        %                     toc
                    end
                    
                    depth = height - com(3);
                    [min_depthTS,min_ind] = min(depth);
                    depths = [min_depthSS,min_depthDS,min_depthTS];
                    [min_depth,type] = min(depths);
                    %                 outputfcn1(type,min_ind);
                    tested = 0;
                    while ~tested && type == 4
                        % it is possible that the TS height based on a stance unconnected to the BG, but only
                        % worth checking when it is dominant
                        [~,L,~] = TSNormals(comsolA(min_ind,:).',p1,p2,p3,nsolA(min_ind,:).');
                        if all(L>0)
                            tested = 1;
                        else
                            printSurf = 0;
                            % The stance found is unfeasible, and therefore
                            % disconnected from the BG. Its region of the
                            % cliff curve is removed, then the next lowest
                            % stance is identified.
                            %                             minNonTSHeight = min([min_depthSS,min_depthDS]);
                            if printSurf
                            CylinderTS2(r,hcm,P,t1,phi,100,1);
                            drawCliff3D(s3_range,phi3_range,height,'k');
                            plot3(0,norm(p3),com(3),'g+','markerSize',6,'linewidth',4)
                            end
                            [phi3_range,s3_range,height,comsolA,nsolA] = cutSectionOffCliffCurve(phi3_range,s3_range,height,comsolA,nsolA,min_ind,[]);
                            if ~isempty(height)
                                [min_depthTS,min_ind] = min(height - com(3));
                            else
                                min_depthTS = inf;
                                tested = 1;
                            end
                            depths = [min_depthSS,min_depthDS,min_depthTS];
                            [min_depth,type] = min(depths);
                        end
                    end
                    
                    if min_depth>0
                        %                 %                     escape_stance_phi3(j,k) = phi3_range(height==min(height));
                        %                 figure(f(i,j))
                        %                 plot(phi3_range,height,'DisplayName',['\theta_1=' num2str(round(t1,2)) ',\phi_1=' num2str(round(phi,2))])
                        escape_height_array(j,k) = min_depth;
                        %                     escape_height_arrayTSA(j,k) = min_heightTSA;
                        escape_height_arrayTS(j,k) = min_depthTS;
                        escape_height_arrayDSRec(j,k) = min_depthDS(1);
                        escape_height_arrayDSCirc(j,k) = min_depthDS(2);
                        escape_height_arraySS(j,k) = min_depthSS;
                        escape_type(j,k) = type;
                    else
                        disp('negative height - computation error')
                    end
                    if(min_depth>maxi) %search for optimally secure grasp
                        max_ds = ds;
                        max_phi = phi;
                        max_theta1 = t1;
                        maxi = min_depth;
                    end
                    
                    %             else
                    %                 fig = Draw3DGrasp(r,hcm*2,hcm,t1,p1,p2,p3);
                    %                 disp('unstable grasp');
                    end
                else
                    disp('unstable grasp');
                end
            end
        end
    end
    %     xlabel('\phi_3');
    %     ylabel('COM height');
    %     title(['h_{cm}=' num2str(hcm)]);
    
    
    %%  approximation error for triple-support approximation
    %   error = (escape_height_arrayTS.'-escape_height_arrayTSA.');
    % rel_error = error./escape_height_arrayTSA.';
    % figure;
    % [~,cf] = contourf(t1_range(1:6),phi_range,error(:,1:6),30);
    % set(cf,'LineColor','none')
    % cb = colorbar;
    % xlabel('$\theta_1$','interpreter','latex','FontSize',20)
    % ylabel('$\phi_1$','interpreter','latex','FontSize',20)
    % set(get(cb,'label'),'string','$\Delta\:U_{3c}$','interpreter','latex','FontSize',20,'rotation',0);
    % title('Approximation error','interpreter','latex','FontSize',20)
    % xlim([min(t1_range),max(t1_range(escape_height_arraySS(:,1)>0))])
    %%
    
    % max_error = max(error,2);
    % max_rel_error = max(rel_error,2);
    %%
  t_total = toc
    %% plot and save figures of depth
    theta = 1;
    printlist = [1,2];
    for c=1:theta
        x_range = t1_range;
%         y_range = phi_range;
        y_range = dsv;
        subFolderName = [pwd '\NewData\Depth_r=' sprintf('%.0e',r) '_ds=' sprintf('%.0e',ds)];
        mkdir(char(subFolderName))     
        if ~isempty(find(printlist==1,1))
            save = 1;
            z_values = escape_height_array;
            z_values(escape_type==0)=-inf;
%             strings = {'$\theta$','$\phi$','$U(q)$',['$h_{12} =$ ' num2str(ds)],'U_with_types'};
            strings = {'$\theta$','$h_{12}$','$U(q)$',['$\phi =$ ' num2str(phi)],'U_with_types'};
            fig = GraspSpaceFig2D(x_range,y_range,z_values,strings,save,subFolderName,escape_type);
        end
        if ~isempty(find(printlist==2,1))
            save = 1;
            z_values = escape_type;
%             strings = {'\theta','\phi','types',['h_{12} = ' num2str(ds)],'U'};
            strings = {'$\theta$','$h_{12}$','$U(q)$',['$\phi =$ ' num2str(phi)],'U_with_types'};
            fig = GraspSpaceFigTypes(x_range,y_range,z_values,strings,save,subFolderName);
        end
        if ~isempty(find(printlist==3,1))
            save = 1;
            error = (escape_height_arrayTS.'-escape_height_arrayTSA.');
            rel_error = error./escape_height_arrayTSA.';
            z_values = rel_error;
            strings = {'\theta','\phi','E_{app}',['Approx. error, h_{12} = ' num2str(ds)],'approx_error'};
            fig = GraspSpaceFig2D(x_range,y_range,z_values,strings,save,subFolderName);
        end
        %%
    end
% end
%%

function [phi3_range,s3_range,height,comsolA,nsolA] = cutSectionOffCliffCurve(phi3_range,s3_range,height,comsolA,nsolA,min_ind,minNonTSHeight)
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
    cutoff2 = find(numer_deriv(min_ind:end)<0,1,'first')+min_ind-1;
    end
    if isempty(cutoff2)
        cutoff2 = length(height);
    end
end
phi3_range(cutoff1:cutoff2) = [];
s3_range(cutoff1:cutoff2) = [];
nsolA(cutoff1:cutoff2,:) = [];
height(cutoff1:cutoff2) = [];
comsolA(cutoff1:cutoff2,:) = [];
end

function outputfcn1(type,min_ind)
switch type
    case 1
        disp('SS');
    case 2
        disp('DS_rec');
    case 3
        disp('DS_circ');
    case 4
        disp('TS');
        disp(min_ind);
end
end

function [p1,p2,com] = grasp_details(r,hcm,phi,t1,ds)
ph1 = phi; ph2 = -phi;
vx = [1;0;0];
vy = [0;cos(t1);sin(t1)];
vz = [0;-sin(t1);cos(t1)];
Rt = [vx vy vz];
com_0 = [0;0;hcm];
com = Rt*com_0;

% position of the contact points
po1 = [0;0;hcm+ds]+r*[sin(ph1);cos(ph1);0];
p1 = Rt*po1;
po2 = [0;0;hcm+ds]+r*[sin(ph2);cos(ph2);0];
p2 = Rt*po2;
end

function fig = GraspSpaceFig2D(x_range,y_range,z_values,strings,save,subFolderName,varargin)
fig = figure;
hold on
% ax(1)=axes;
[~,cf] = contourf(x_range,y_range,z_values.',30);
ax(1) = gca;
% set(cf,'LineColor','none')
cb = colorbar;
x_string = strings{1};
y_string = strings{2};
z_string = strings{3};
xlabel(x_string,'interpreter','latex','FontSize',20)
ylabel(y_string,'interpreter','latex','FontSize',20)
set(get(cb,'label'),'string',z_string,'interpreter','latex','FontSize',20,'rotation',0);

if nargin>=7
    escape_type = varargin{1};
    ax(2) = axes;
    contour(x_range,y_range,escape_type.',[1.5,2.5,3.5],'-k','linewidth',4);
    % colormap(ax(2),[0 0 0])
    pos = get(ax(1),'Position');
    set(ax(2),'Position',pos);
    get(ax(2),'Position')
    set(ax(2),'Color','none',...
        'XTick',[],...
        'YTick',[]);
    
    escape_type_v = reshape(escape_type,1,[]).';
    Nx = length(x_range);
%     for i=1:4
%         relevant_escape = find(escape_type_v==i);
%         if ~isempty(relevant_escape)
% %             xmean = mean(x_range(floor((relevant_escape-1)/Nx+1)));
% %             ymean = mean(y_range(mod(relevant_escape,Nx)+1));
% %             text(xmean,ymean,num2str(i),'fontsize',24);
%         end
%     end
%     xlim(ax(2),[min(x_range),max(x_range(z_values(:,1)>0))])
end
hold off
% xlim(ax(1),[min(x_range),max(x_range(z_values(:,1)>0))])


title_string = strings{4};
title(title_string,'interpreter','latex','fontsize',20)
if save == 1
    filename = strings{5};
savefig(fig,[subFolderName filename])
saveas(fig,[subFolderName filename '.eps'],'epsc')
saveas(fig,[subFolderName filename '.png'],'png')
end
end

function fig = GraspSpaceFigTypes(x_range,y_range,z_values,strings,save,subFolderName,varargin)
fig = figure;
hold on
% ax(1)=axes;
[~,cf] = contourf(x_range,y_range,z_values.',[1,2,3,4]);
ax(1) = gca;
% set(cf,'LineColor','none')
cb = colorbar;
x_string = strings{1};
y_string = strings{2};
z_string = strings{3};
xlabel(x_string,'interpreter','latex','FontSize',20)
ylabel(y_string,'interpreter','latex','FontSize',20)
set(get(cb,'label'),'string',z_string,'interpreter','latex','FontSize',20,'rotation',0);

hold off
% xlim(ax(1),[min(x_range),max(x_range(z_values(:,1)>0))])


title_string = strings{4};
title(title_string,'interpreter','latex','fontsize',20)
if save == 1
    filename = strings{5};
savefig(fig,[subFolderName filename])
saveas(fig,[subFolderName filename '.eps'],'epsc')
saveas(fig,[subFolderName filename '.png'],'png')
end
end