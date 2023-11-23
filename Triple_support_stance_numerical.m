function [height,comsol,angles,nsol] = Triple_support_stance_numerical(r,hcm,t,P,t1_feq1,phi3,d3,options) %#ok<INUSL>
comsol = [];angles = [];nsol = [];
h = cos(t); %height of com at basket grasp
% t1_feq1 = @(x) t1_feq2(x,d3);
height = inf;
b = (P(3,1)- P(3,3));
a = -(P(2,1) - P(2,3));
mag = sqrt(a^2+b^2);
numer = abs(d3*sin(phi3));
epsi = 1E-3;
%there is a range of t1's that will lead to complex values, exclude it
if mag<numer
    %no solution exists (formula for t2 returns only complex values)
    return
end
phase = atan2(b,a);
angle = abs(asin(d3*sin(phi3)/mag));

bounds = [-phase+angle+epsi,pi()-phase-angle-epsi];

steps = 100; guesses = initial_guess(t1_feq1,bounds,steps);
bounds = [pi()-phase+angle+epsi,2*pi()-phase-angle-epsi];
steps = 100; guesses = [guesses;initial_guess(t1_feq1,bounds,steps)];

if isempty(guesses)
    return
end

p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
% p1x = p1(1);
p1y = p1(2);
p1z = p1(3);
p3y = p3(2);
p3z = p3(3);

%     tic
for j=1:size(guesses,1)
    guess = guesses(j,:);
    try
        [t1_sol,~,~,~] = fzero(t1_feq1,guess,options);
        t1_sol = wrapToPi(t1_sol);
%         if abs(t1_sol)>pi()/2
%             continue;
%         end
        
    catch me %#ok<NASGU>
        disp('error in fzero')
    end
    
    
    
    t2_sol =  -asin((d3*sin(phi3))./(p1z*cos(t1_sol) - p3z*cos(t1_sol) - p1y*sin(t1_sol) + p3y*sin(t1_sol)));
    t1t = t1_sol;
    t2t = t2_sol;
    Rxsol = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
    Rysol = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
    Rtsol = Rxsol*Rysol;
    nsol = Rtsol*[0;0;1];
    centersol = p3 + Rxsol*Rysol*[d3*sin(phi3);-d3*cos(phi3);0];
    comsol = centersol+nsol*hcm;
    height = comsol(3);
    
%          fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,t1t,t2t,phi3,p1,p2,p3);
    %% check if center of mass outside the triangle defined by the fingers
    if 0
        xq = comsol(1); %#ok<UNRCH>
        yq = comsol(2);
        xv = [p1(1) p2(1) p3(1)];
        yv = [p1(2) p2(2) p3(2)];
        if ~inpolygon(xq,yq,xv,yv)
            [~,L,~] = TSNormals(r,hcm,comsol,phi3,p1,p2,p3,[t1t,t2t]);
            if ~all(L>0)
                height = inf;
            else
                disp('feasible eq in expected infeasible eq region');
            end % 
        end
    end
    
    %% check if the side fingers crossed the bases edge
    dist1 = (p1-centersol).'*nsol;
    dist2 = (p2-centersol).'*nsol;
    if any([dist1,dist2]<0)||any([dist1,dist2]>2*hcm)
        height=inf;
        %         break_cond = 1;
    end
    %     if any([dist1,dist2]>2*hcm)
    %         disp('fingers too high')
    %     end
    %% check if axis in solution "crossed over" point between fingers
    % if it did, since the initial guess is supposed to be closer to the correct solution, that implies that the relevant solution doesn't exist.
    p12 = (p1+p2)/2;
    comsolp = comsol;
    comsolp(1) = 0;
    centersolp = centersol;
    centersolp(1) = 0;
    nsolp = nsol;
    nsolp(1) = 0;
    s = (comsolp-centersolp).'*(centersolp-p12)/norm(comsolp-centersolp)^2;
    arm = (centersolp-s*(comsolp-centersolp))-p12;
    torque = cross(arm,nsolp);
    if torque(1)<0
        height = inf;
        %         [dist1,dist2]
    end
    
    %   fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
%         if h-height>1E-4
%              disp(['possible instability, negative depth = '  num2str(h-height)])
%         end
    
    if height ~= inf
        angles = [t1t,t2t];
        return
    end
end
end

function debug_plotter(d3,phi3,P) %#ok<DEFNU>
p1 = P(:,1);
% p2 = P(:,2);
p3 = P(:,3);
% p1x = p1(1);
p1y = p1(2);
p1z = p1(3);
p3y = p3(2);
p3z = p3(3);
t1_range = linspace(-pi(),pi(),10000);
t2_range = -asin((d3*sin(phi3))./((p1z- p3z)*cos(t1_range) - (p1y - p3y)*sin(t1_range)));
figure; plot(t1_range,(d3*sin(phi3))./((p1z- p3z)*cos(t1_range) - (p1y - p3y)*sin(t1_range)))

figure; plot(t1_range,imag(t2_range));
% figure; plot(t1_range,real(t2_range));
end

function guesses = initial_guess(t1_feq11,bounds,steps)
plot_t = linspace(bounds(1),bounds(2),steps);val = t1_feq11(plot_t);
dif_val = diff(val);
dif_val = [dif_val(1) dif_val];
guesses = zeros(3,2);
c = 1;
while length(val)>2
    try
        end_ind = find(dif_val*dif_val(1)<=0,1,'first')-1;
        if isempty(end_ind)
            end_ind = length(dif_val);
        end
        if end_ind==0
            end_ind = 1;
        end
        if val(1)*val(end_ind)<0
            guesses(c,:) = plot_t([1,end_ind]);
            c = c+1;
        else
            %if endpoint is close to zero, worth it to check if a finer
            %discretization between the previous and the next point will find a
            %zero crossing
            if abs(val(end_ind))<abs(val(1))
                if end_ind==length(val)
                    end_ind = end_ind-1;
                end
                plot_t_zoom = linspace(plot_t(end_ind-1),plot_t(end_ind+1),steps);val_zoom = t1_feq11(plot_t_zoom);
                neg_ind = find(val_zoom*val(1)<0,1,'first');
                if ~isempty(neg_ind)
                    guesses(c,:) = [plot_t(1),plot_t_zoom(neg_ind)];
                    c = c+1;
                end
            end
        end
        plot_t(1:end_ind) = [];
        val(1:end_ind) = [];
        dif_val(1:end_ind) = [];
    catch me %#ok<NASGU>
        disp('problem with guess range finder')
    end
end
guesses(guesses*[1;1]==0,:) = [];
end





























%% old versions of functions 
% function [guess1,guess2] = initial_guess_old(t1_feq11,bounds,steps)
%     
% %     val1 = inf*ones(1,steps);
%     plot_t = zeros(1,steps);
% %     bounds = [-pi(),pi()];
%     if bounds(1)>bounds(2)
% %         disp('sus')
%         guess1 = [];
%         guess2 = [];
%         return
%     end
%     peaks = [-inf,-inf];
%     while any(peaks<0)
%         plot_t = linspace(bounds(1),bounds(2),steps);val = t1_feq11(plot_t);%figure; plot(plot_t,val1);
%         if val(1)<0
%             next_pos = find(val>0,1,'first');
%             next_neg = find(val(next_pos:end)<0,1,'first');
%             if ~isempty(next_neg)
%                 val(1:next_pos) = [];
%                 plot_t(1:next_pos) = [];
%             end
%         end
%         
%         if any(imag(val)~=0)
%             disp('exclusion of complex angles failed')
%         end
% 
%         if any(val<0)
%             if val(1)<0
%                 ind_guess = find(val>0,1,'first');
%                 break
%             else
%                 ind_guess = find(val<0,1,'first');
%                 break
%             end
%         end
% [~, inds] = findpeaks(-real(val));
%         if length(inds)>1
%             ind_guess = inds(abs(inds-0.5*length(val))==min(abs(inds-0.5*length(val))));
%         else
%             ind_guess = inds; 
%             if isempty(inds)
%                 ind_guess = find(val<0,1,'first');
%             end 
%         end
%         if ~isempty(ind_guess)
%         peaks = -val(ind_guess);
%         prev_bounds = bounds;
%         bounds = [plot_t(ind_guess-1),plot_t(ind_guess+1)];
%         else
%             guess1 = [];
%             guess2 = [];
%             return
%         end
%     end
% %         bounds_dif = diff(prev_bounds);
% %         bounds_dif = min(bounds_dif,pi()/4);
%         theta_neg = real(plot_t(ind_guess));
%         if val(1)>0
%             %         theta_neg = theta_neg + pi()*(abs(theta_neg)>pi()/2);
%             %         guess1 = [theta_neg-bounds_dif/2,theta_neg];
%             %         guess2 = [theta_neg,theta_neg+bounds_dif/2];
%             guess1_st = plot_t(find(val(1:ind_guess)>0,1,'last'));
%             guess1 = [guess1_st,theta_neg];
%             guess2_end = plot_t(ind_guess -1 + find(val(ind_guess:end)>0,1,'first'));
%             guess2 = [theta_neg,guess2_end];
%         else
%             guess1_st = plot_t(find(val(1:ind_guess)<0,1,'last'));
%             guess1 = [guess1_st,theta_neg];
%             guess2_end = plot_t(ind_guess -1 + find(val(ind_guess:end)<0,1,'first'));
%             guess2 = [theta_neg,guess2_end];
%         end
%         
% end

