function varargout = CylinderTSCurve3(r,d3,hcm,P,t,phi,N,prev_sol,feq2,t2_sol2)

% phi3_range = linspace(-pi(),0,N);%linspace(-1.8,-0.8,N);%


% options = optimoptions('fsolve','OptimalityTolerance',1E-8,'Display','none','FunctionTolerance',1E-8);
options0 = optimset('Display','none','TolX',1e-8);

% [t1s1,t1s2] = TS_symmetrical2(r,d3,P,phi); %t1 for phi3=0,pi
phi3_range = linspace(-pi(),0,N);
[circle_height,circleTheta1,circleTheta2,centersol,comsol,nsol] = computeHalfCircle(phi3_range,r,d3,hcm,P,t,N,options0,prev_sol,feq2,t2_sol2);


%%

switch nargout
    case 0
        % save('circle4.mat','phi3_range','circle_height','circleTheta1','circleTheta2');
        figure
        hold on
        plot(phi3_range,circle_height)
        % plot(circleTheta1,circle_height)
        xlabel('\phi_3')
        ylabel('COM height');
        fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,minCOM,minTheta1,minTheta2,p1,p2,p3);
    case 1
        varargout{1} = circle_height;
    case 2
        varargout{1} = phi3_range;
        varargout{2} = circle_height;
    case 4
        varargout{1} = phi3_range;
        varargout{2} = circle_height;
        varargout{3} = circleTheta1;
        varargout{4} = circleTheta2;
    case 3
        varargout{1} = minCOM;
        varargout{2} = minTheta1;
        varargout{3} = minTheta2;
end
end

function [circleHeight,circleTheta1,circleTheta2,centersolA,comsolA,nsolA] = computeHalfCircle(phi3_range,r,d3,hcm,P,t,N,options,prev_sol,feq2,t2_sol2)
circleHeight = inf*ones(size(phi3_range));
circleTheta1 = zeros(size(phi3_range));
circleTheta2 = inf*ones(size(phi3_range));
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
minh = inf;
s3 = norm(p3);
epsi = 1E-2;
h = CylinderComPos(r,hcm,P,t,0,s3);
h= cos(t);
% supportX = [p1(1), p2(1), p3(1)];
% supportY = [p1(2), p2(2), p3(2)];
centersolA = zeros(N,3);
comsolA = centersolA;
nsolA = centersolA;
t1_sol = [];
break_cond = [];
for i=1:length(phi3_range)
    phi3 = phi3_range(i);
    t1_feq21 = feq2{1};
    t1_feq22 = feq2{2};
    t1_feq11 = @(x) t1_feq21(x,phi3);
    [height,comsol,angles,nsol] = Triple_support_stance_numerical(r,hcm,t,P,t1_feq11,phi3,d3,options);
    circleHeight(i) = height;
    if ~isempty(angles)
    circleTheta1(i) = angles(1);
    circleTheta2(i) = angles(2);
    centersolA(i,:) = (comsol-nsol).';
    comsolA(i,:) = comsol.';
    nsolA(i,:) = nsol.';
    end
    for non=[]
%     t1_feq12 = @(x) t1_feq22(x,phi3);
%     feq1 = @(x) feq2(x,phi3);
    t2_sol1 = @(x) t2_sol2(x,phi3);
    
%     x = sym('x',[2,1]);
%     t1t = x(1);
%     t2t = x(2);
    
%     p2t = p2-p3;
%     p1t = p1-p3;
%     Rx = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
%     Ry = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
%     Rt = Rx*Ry;
%     nt = [0;0;1];
%     Rz = [cos(phi3) -sin(phi3) 0;sin(phi3) cos(phi3) 0;0 0 1];
%     center = Rt*Rz*[0 -r 0].';
%     n = Rt*nt;
%     p1p = p1t-(p1t.'*n)*n;
%     p2p = p2t-(p2t.'*n)*n;
    % rad1 = p1p-center;
%     rad1 = cross(p1t-center,n);
    % rad2 = p2p-center;
%     rad2 = cross(p2t-center,n);
%     eq1 = rad1.'*rad1-r^2;
%     eq2 = rad2.'*rad2-r^2;
%     eq = [eq1;eq2];
%     eq = simplify(eq)
%     eqex = expand(eq);
%     eq3 = simplify(eq(1)-eq(2));
    % t1_sol = solve(eq3,x(1))
%     t2_sol = solve(eq3,x(2), 'Real', true)
%     if i==1
%     guess = prev_sol(:,i); %take as initial guess the solution from the same angle, one step close to the origin
%     else
%         guess = t1_sol;
%         if isempty(t1_sol)
%             guess = [t;t];
%         end
%         
%     end
%     guess(guess==inf)=t;
%     guess(imag(guess)~=0)=t;
%     guess = [-1;-1];
% tic



    bounds = [-pi()/2 pi()/2];
    b = (P(3,1)- P(3,3));
    a = -(P(2,1) - P(2,3));
    mag = sqrt(a^2+b^2);
    numer = abs(d3*sin(phi3));
    epsi = 1E-3;
 %there is a range of t1's that will lead to complex values, exclude it
 if mag<numer
     %no solution exists (formula for t2 returns only complex values)
     continue
 end
 phase = atan2(b,a);
angle = abs(asin(d3*sin(phi3)/mag));

bounds = [-phase+angle+epsi,pi()-phase-angle-epsi];

steps = 100; guesses = initial_guess(t1_feq1,bounds,steps);
bounds = [pi()-phase+angle+epsi,2*pi()-phase-angle-epsi];
steps = 100; guesses = [guesses;initial_guess(t1_feq1,bounds,steps)];
 
 if 0
        phase = atan2(b,a);
        angle = abs(asin(d3*sin(phi3)/mag));
        bounds = [wrapToPi(-phase+angle+epsi),min(pi()/2,wrapToPi(pi()-phase-angle-epsi))];
        steps = 100;[guess1,guess2] = initial_guess(t1_feq11,bounds,steps);
        if isempty([guess1,guess2])
        bounds = [max(-pi()/2,wrapToPi(-phase-angle+epsi)),wrapToPi(-phase-angle-epsi)];
        steps = 100;[guess1,guess2] = initial_guess(t1_feq11,bounds,steps);
        end
 end
%     else
%         steps = 100;[guess1,guess2] = initial_guess(t1_feq11,bounds,steps);
%     end
% toc_ig = toc
if isempty([guess1,guess2])
    continue
end
try
%     tic
    [t1_sol1,value,flag,~] = fzero(t1_feq11,guess1,options);
    [t1_sol2,value,flag,~] = fzero(t1_feq11,guess2,options);
    t1_sol = [t1_sol1;t1_sol2];
    t1_sol = wrapToPi(t1_sol);
    
    if t1_sol1<-pi()/2 ||t1_sol1>pi()/2 || any(imag(t1_sol)~=0)
        disp('maybe need to repeat fzero with diff. initial guess');
    end

catch me
         disp('error in fzero')
end
    
%     catch
%         disp('error in fsolve')
%     end
    if flag < 1 || norm(value)>1E-2
        t1_sol = [];
        continue
    end
    t2_sol = feval(t2_sol1,t1_sol);
 for j=1:2   
    t1t = t1_sol(j);
    t2t = t2_sol(1);
    Rxsol = [1 0 0; 0 cos(t1t) -sin(t1t); 0 sin(t1t) cos(t1t)];
    Rysol = [cos(t2t) 0 sin(t2t);0 1 0; -sin(t2t) 0 cos(t2t)];
    Rtsol = Rxsol*Rysol;
    nsol = Rtsol*[0;0;1];
    nsolA(i,:) = nsol.';
    centersol = p3 + Rxsol*Rysol*[d3*sin(phi3);-d3*cos(phi3);0];
    centersolA(i,:) = centersol.';
    comsol = centersol+nsol*hcm;
    comsolA(i,:) = comsol.';
    ht = comsol(3);
    circleHeight(i) = ht;
    circleTheta1(i) = t1t;
    circleTheta2(i) = t2t;
    %         fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,t1t,t2t,phi3,p1,p2,p3);
    %% check if center of mass outside the triangle defined by the fingers
    if 0
    xq = comsol(1);
    yq = comsol(2);
    xv = [p1(1) p2(1) p3(1)];
    yv = [p1(2) p2(2) p3(2)];
    if ~inpolygon(xq,yq,xv,yv)
        
        [~,L,~] = TSNormals(r,hcm,comsol,phi3,p1,p2,p3,[t1t,t2t]);
        if all(L>0)
%             disp('feasible eq in expected infeasible eq region');
        else
            circleHeight(i) = inf;
        end
    end
    end
    
    %% check if the side fingers crossed the bases edge
    dist1 = (p1-centersol).'*nsol;
    dist2 = (p2-centersol).'*nsol;
    if any([dist1,dist2]<0)||any([dist1,dist2]>2*hcm)
        circleHeight(i)=inf;
        break_cond = 1;
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
        circleHeight(i) = inf;
%         [dist1,dist2]
    end
    x = d3*sin(phi3);
    y = d3*cos(phi3);

%   fig = Draw3DGrasp(r,hcm*2,hcm,t,p1,p2,p3);Draw3DEscape(r,hcm*2,hcm,comsol,sol(1),sol(2),phi3,p1,p2,p3);
%     if h-circleHeight(i)>1E-4
%          disp(['possible instability, negative depth = '  num2str(h-circleHeight(i))])
%     end
    
    if circleHeight(i) ~= inf
        break
    end
 end
 if d3 ==0
     circleHeight(1:end) = circleHeight(i);
     break
 end
    end
end
end


function debugger1(d3,phi3,P)
p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p1x = p1(1);
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


function [guess1,guess2] = initial_guess_old(t1_feq11,bounds,steps)
    
%     val1 = inf*ones(1,steps);
    plot_t = zeros(1,steps);
%     bounds = [-pi(),pi()];
    if bounds(1)>bounds(2)
%         disp('sus')
        guess1 = [];
        guess2 = [];
        return
    end
    peaks = [-inf,-inf];
    while any(peaks<0)
        plot_t = linspace(bounds(1),bounds(2),steps);val1 = t1_feq11(plot_t);%figure; plot(plot_t,val1);
        if val1(1)<0
            next_pos = find(val1>0,1,'first');
            next_neg = find(val1(next_pos:end)<0,1,'first');
            if ~isempty(next_neg)
                val1(1:next_pos) = [];
                plot_t(1:next_pos) = [];
            end
        end
        
        if any(imag(val1)~=0)
            disp('exclusion of complex angles failed')
        end

        if any(val1<0)
            if val1(1)<0
                ind_guess = find(val1>0,1,'first');
                break
            else
                ind_guess = find(val1<0,1,'first');
                break
            end
        end
[~, inds] = findpeaks(-real(val1));
        if length(inds)>1
            ind_guess = inds(abs(inds-0.5*length(val1))==min(abs(inds-0.5*length(val1))));
        else
            ind_guess = inds; 
            if isempty(inds)
                ind_guess = find(val1<0,1,'first');
            end 
        end
        if ~isempty(ind_guess)
        peaks = -val1(ind_guess);
        prev_bounds = bounds;
        bounds = [plot_t(ind_guess-1),plot_t(ind_guess+1)];
        else
            guess1 = [];
            guess2 = [];
            return
        end
    end
%         bounds_dif = diff(prev_bounds);
%         bounds_dif = min(bounds_dif,pi()/4);
        theta_neg = real(plot_t(ind_guess));
        if val1(1)>0
            %         theta_neg = theta_neg + pi()*(abs(theta_neg)>pi()/2);
            %         guess1 = [theta_neg-bounds_dif/2,theta_neg];
            %         guess2 = [theta_neg,theta_neg+bounds_dif/2];
            guess1_st = plot_t(find(val1(1:ind_guess)>0,1,'last'));
            guess1 = [guess1_st,theta_neg];
            guess2_end = plot_t(ind_guess -1 + find(val1(ind_guess:end)>0,1,'first'));
            guess2 = [theta_neg,guess2_end];
        else
            guess1_st = plot_t(find(val1(1:ind_guess)<0,1,'last'));
            guess1 = [guess1_st,theta_neg];
            guess2_end = plot_t(ind_guess -1 + find(val1(ind_guess:end)<0,1,'first'));
            guess2 = [theta_neg,guess2_end];
        end
        
end

