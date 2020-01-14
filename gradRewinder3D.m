function [grad_waveform,k_traj,slewRate] = gradRewinder3D(g0,g_des,k0,k_des,GMAX,SMAX,TS,nucl)
%
%  Attemtps to link two points in the 3D k-space in the shortest possible 
%  time with prescribed initial and final gradients and constraints on the 
%  maximal gradient amplitude and slew-rate. 
%
%  ************************************************************************
%  !!!Not working perfectly, with some parameters gives erroneous result!!!
%  ************************************************************************
%
%   INPUT
%       g0    = current gradient position [G/cm]
%       g_des = current gradient position [G/cm]
%       k0    = current k-space position
%       GMAX  = maximal gradient amplitude [G/cm]
%       SMAX  = maximal slew-rate [G/cm/s]
%       TS    = sampling time [s]
%
%   OUTPUT
%       grad_waveform = rewinder gradient trajectory [G/cm]
%       k_traj        = rewinder k-space trajectory [1/m]
%       slewRate      = slew-rate [G/cm/s]   
%
%   Vencel Somai -> vs460@cam.ac.uk

if strcmp(nucl,'1H') 
    % 1H gyromagnetic-ratio
    gamma = 4258;
elseif strcmp(nucl,'13C')
    % 13C gyromagnetic-ratio
    gamma = 1071;
elseif strcmp(nucl,'2H')
    % 2H gyromagnetic-ratio
    gamma = 653.6;
else
    msg = 'nucl should be 1H or 13C\n';
    error(msg);
end

% then there are two cases: slew-rate limited and amplitude limited
% first: slew-rate limited
delta_k = k_des-k0;
abs_delta_k = abs(delta_k);
signsOfSlewRates = ((abs(g0)+abs(g_des))/2*TS < abs_delta_k);
if sum(abs((GMAX-abs(g0))/SMAX.*(g0+sign(delta_k)*GMAX)/2 + (GMAX-abs(g_des))/SMAX.*(g_des+sign(delta_k)*GMAX)/2) > abs(delta_k)) == 3;

    % calculate the maximal gradient amplitude needed so that the triangle
    % gradient has the given area = k-space location
    g_max = sign(delta_k).*abs(sqrt((signsOfSlewRates.*delta_k/gamma/100*SMAX+(g0.^2+g_des.^2)/2)));
    num_of_t_pts_up = ceil(abs(g_max-g0)/SMAX/TS);
    num_of_t_pts_down = ceil(abs(g_max-g_des)/SMAX/TS);
    % calculate the number of points needed to remain below the maximal
    % slew-rate
    g_max = (2*delta_k/gamma/TS/100 - g0.*(num_of_t_pts_up+1) - g_des.*num_of_t_pts_down)./(num_of_t_pts_up+1+num_of_t_pts_down);

    grad_waveform = zeros(max(num_of_t_pts_up+num_of_t_pts_down+1),3);
    % rising part x-axis
    grad_waveform(1:num_of_t_pts_up(1)+1,1) =...
            [linspace(g0(1),g_max(1),num_of_t_pts_up(1)+1)'];
    % decaying part x-axis
    grad_waveform(num_of_t_pts_up(1)+2:num_of_t_pts_up(1)+num_of_t_pts_down(1)+1,1)  =...
            [linspace(g_max(1),g_des(1),num_of_t_pts_down(1))'];
    % rising part y-axis
    grad_waveform(1:num_of_t_pts_up(2)+1,2) =...
             [linspace(g0(2),g_max(2),num_of_t_pts_up(2)+1)'];

    % decaying part y-axis
    grad_waveform(num_of_t_pts_up(2)+2:num_of_t_pts_up(2)+num_of_t_pts_down(2)+1,2)  =...
            [linspace(g_max(2),g_des(2),num_of_t_pts_down(2))'];
    % rising part z-axis
    grad_waveform(1:num_of_t_pts_up(3)+1,3) =...
            [linspace(g0(3),g_max(3),num_of_t_pts_up(3)+1)'];
    % decaying part z-axis
    grad_waveform(num_of_t_pts_up(3)+2:num_of_t_pts_up(3)+num_of_t_pts_down(3)+1,3)  =...
            [linspace(g_max(3),g_des(3),num_of_t_pts_down(3))'];
    
    idx_recalc = find(abs(sum(grad_waveform,1)*TS*gamma*100-delta_k)>0.01*abs(delta_k));
    g_max(idx_recalc) = sign(delta_k(idx_recalc)).*abs(sqrt((-signsOfSlewRates(idx_recalc).*delta_k(idx_recalc)*2/gamma/100*SMAX+(g0(idx_recalc).^2+g_des(idx_recalc).^2))/2));

    % calculate the number of points needed to remain below the maximal
    % slew-rate
    num_of_t_pts_up = ceil(abs(g_max-g0)/SMAX/TS);
    num_of_t_pts_up = max(num_of_t_pts_up);
    num_of_t_pts_down = ceil(abs(g_max-g_des)/SMAX/TS);
    num_of_t_pts_down = max(num_of_t_pts_down);
    % recalculate g_max so that all the axes are rewinded at the same time
    slewRateMet = false;
    while ~slewRateMet
        g_max = (2*delta_k/gamma/TS/100 - g0*(num_of_t_pts_up+1) - g_des*num_of_t_pts_down)/(num_of_t_pts_up+1+num_of_t_pts_down);
        % calculate the triangle gradient waveform
        % rising part
        grad_waveform(1:num_of_t_pts_up+1,:) =...
            [linspace(g0(1),g_max(1),num_of_t_pts_up+1)',linspace(g0(2),g_max(2),num_of_t_pts_up+1)',linspace(g0(3),g_max(3),num_of_t_pts_up+1)'];
        % decaying part
        grad_waveform(num_of_t_pts_up+2:num_of_t_pts_up+num_of_t_pts_down+1,:)  =...
            [linspace(g_max(1),g_des(1),num_of_t_pts_down)',linspace(g_max(2),g_des(2),num_of_t_pts_down)',linspace(g_max(3),g_des(3),num_of_t_pts_down)'];
        if max(max(abs(diff(grad_waveform,1)/TS))) > SMAX
            num_of_t_pts_up = num_of_t_pts_up + 1;
            num_of_t_pts_down = num_of_t_pts_down + 1;
        else
            slewRateMet = true;
        end
    end
    
% second: amplitude limited case
else
    [maxValue,axis] = max(abs(delta_k/gamma/100 - (GMAX-abs(g0))/SMAX.*(g0+sign(delta_k)*GMAX)/2 - (GMAX-abs(g_des))/SMAX.*(g_des+sign(delta_k)*GMAX)/2));
    % how many time points are needed to reach maximal gradient amplitude
    num_of_rise_pts = ceil(abs(sign(delta_k(axis))*GMAX-(g0(axis)))/SMAX/TS);
    % how many time points are needed to reach the desired gradient
    % amplitude from the maximal
    num_of_decay_pts = ceil(abs(sign(delta_k(axis))*GMAX-g_des(axis))/SMAX/TS);
    slewRateMet = false;
    while ~slewRateMet
        % the area covered by the rising part
        k_area = num_of_rise_pts.*(sign(delta_k(axis))*GMAX+g0(axis))*TS*gamma*100/2;
        % the area covered by decaying part
        k_area = k_area + num_of_decay_pts.*(sign(delta_k(axis))*GMAX+g_des(axis))*TS*gamma*100/2;
        % remaining area where the gradient amp = GMAX
        k_remaining = delta_k(axis) - k_area;
        % calculating how many points it requires
        num_of_flat_pts = ceil(abs(k_remaining/GMAX/TS/gamma/100));
        % calculating the total number of time points
        num_of_t_pts = num_of_rise_pts + num_of_flat_pts + num_of_decay_pts;
        % calculating the maximal gradient which is somewhat smaller than GMAX
        % due to rounding
        g_max = (delta_k(axis)/gamma/100 - g0(axis)/2*num_of_rise_pts*TS - g_des(axis)/2*num_of_decay_pts*TS)/(num_of_rise_pts*TS/2 + num_of_decay_pts*TS/2 + num_of_flat_pts*TS);
        
        % calculating axis - gradient waveform
        grad_waveform(1:num_of_rise_pts,axis) = linspace(g0(axis),g_max,num_of_rise_pts)';
        grad_waveform(num_of_rise_pts+1:num_of_rise_pts+num_of_flat_pts,axis) = linspace(g_max,g_max,num_of_flat_pts)';
        grad_waveform(num_of_rise_pts+num_of_flat_pts+1:num_of_rise_pts+num_of_flat_pts+num_of_decay_pts,axis) = linspace(g_max,g_des((axis)),num_of_decay_pts)';
        if max(max(abs(diff(grad_waveform(:,axis),1)/TS))) > SMAX
            num_of_rise_pts = num_of_rise_pts + 1;
            num_of_decay_pts = num_of_decay_pts + 1;
        else
            slewRateMet = true;
        end
    end  
    
    % calculating the gradient waveform on the remaining two axes
    remaining_axes = find([1,2,3]~=axis);
    g_max = (delta_k(remaining_axes)/gamma/100 - g0(remaining_axes)/2*num_of_rise_pts*TS - g_des(remaining_axes)/2*num_of_decay_pts*TS)/(num_of_rise_pts*TS/2 + num_of_decay_pts*TS/2 + num_of_flat_pts*TS);
    %g_max = (delta_k(remaining_axes)/gamma/100/TS*SMAX*2+g0(remaining_axes).^2+g_des(remaining_axes).^2)./(num_of_t_pts*SMAX*2+g0(remaining_axes)+g_des(remaining_axes));
    num_of_rise_pts = ceil(abs(g_max-g_des(remaining_axes))/SMAX/TS);
    num_of_decay_pts = ceil(abs(g_max-g0(remaining_axes))/SMAX/TS);
    num_of_flat_pts = num_of_t_pts - num_of_rise_pts - num_of_decay_pts;
    
    slewRateMet = false;
    while ~slewRateMet
        g_max = (delta_k(remaining_axes)/gamma/100 - g0(remaining_axes)/2.*num_of_rise_pts*TS - g_des(remaining_axes)/2.*num_of_decay_pts*TS)./(num_of_rise_pts*TS/2 + num_of_decay_pts*TS/2 + num_of_flat_pts*TS);
        % calculate the triangle gradient waveform
        % rising part
        % corresponding waveform 1
        grad_waveform(1:num_of_rise_pts(1),remaining_axes(1)) = linspace(g0(remaining_axes(1)),g_max(1),num_of_rise_pts(1))';
        grad_waveform(num_of_rise_pts(1)+1:num_of_rise_pts(1)+num_of_flat_pts(1),remaining_axes(1)) = linspace(g_max(1),g_max(1),num_of_flat_pts(1))';
        grad_waveform(num_of_rise_pts(1)+num_of_flat_pts(1)+1:num_of_rise_pts(1)+num_of_flat_pts(1)+num_of_decay_pts(1),remaining_axes(1)) = linspace(g_max(1),g_des(remaining_axes(1)),num_of_decay_pts(1))';
        % corresponding waveform 2
        grad_waveform(1:num_of_rise_pts(2),remaining_axes(2)) = linspace(g0(remaining_axes(2)),g_max(2),num_of_rise_pts(2))';
        grad_waveform(num_of_rise_pts(2)+1:num_of_rise_pts(2)+num_of_flat_pts(2),remaining_axes(2)) = linspace(g_max(2),g_max(2),num_of_flat_pts(2))';
        grad_waveform(num_of_rise_pts(2)+num_of_flat_pts(2)+1:num_of_rise_pts(2)+num_of_flat_pts(2)+num_of_decay_pts(2),remaining_axes(2)) = linspace(g_max(2),g_des(remaining_axes(2)),num_of_decay_pts(2))';
        if max(max(abs(diff(grad_waveform(:,remaining_axes),1)/TS))) > SMAX
            num_of_rise_pts = num_of_rise_pts + 1;
            num_of_decay_pts = num_of_decay_pts + 1;
        else
            slewRateMet = true;
        end
    end
    
    dummyVariable = ones(size(grad_waveform,1));
end

% calculate the corresponding k-space trajectory
k_traj(:,1) = cumsum(grad_waveform(:,1)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
k_traj(:,2) = cumsum(grad_waveform(:,2)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
k_traj(:,3) = cumsum(grad_waveform(:,3)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio


% calculate the slew-rate
slewRate = diff(grad_waveform,1)/TS;
end



