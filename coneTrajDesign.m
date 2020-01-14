function [gw,k_traj,density,N_cones,coneLengths,rewinderBoolean,kv_traj] = coneTrajDesign(RES,FOV,Nint,theta,nucl,MAXLEN,TS,SMAX,GMAX,cvx_exists)

% Cone trajectory design for 13C imaging using the functions from [1]
%
%   [1]: Paul T. Gurney,1* Brian A. Hargreaves,2 and Dwight G. Nishimura1,
%        Design and Analysis of a Practical 3D Cones Trajectory, MRM 2006
% 
% which are available for download at:
%   http://mrsrl.stanford.edu/~ptgurney/cones.php
%
%   Input:
%
%       RES        = isotropic resolution [mm]
%       FOV        = isotropic FOV [cm]
%       Nint       = number of interlieves for the individual cones
%       theta      = range of polar angles
%       MAXLEN     = Maximum length output allowed
%       TS         = Sampling period in seconds
%       SMAX       = Slew rate [G/cm/s]
%       GMAX       = Max Gradient Amplitude [G/cm]
%       cvx_exists = boolean whether cvx package is available or not 
%
%   Output:
%
%        g         = gradient waveform
%        k         = k-space trajectory
%
%   The script generates a single gradient waveform which is then scaled
%   down(!) to each cone, therefore not exceeding hardware limits.
%
%   For rewinding the cones two oprions are available. One is based on
%   symbolic calculations which is fast to compute, is very close to 
%   time-optimal but can give unnecessarily long rewinder in case of some 
%   design parameters. The other is based on the cvx solver (requires cvx 
%   package, see the description in the function) and slower to compute the
%   result but gives time-optimal rewinder for any set of parameters 
%   (at least tested for many sets of parameters and worked fine). 
%   It is worth checking both for a gicen design.
%
%   Vencel Somai -> vs460@cam.ac.uk


%7T Agilent preclinical specifications
if (nargin<6)
    MAXLEN = 10000;
end
if (nargin<7)
    TS = 0.000004;
end
if (nargin<8)
    SMAX = 20000;
end
if (nargin<9)
    GMAX = 4;
end
if (nargin<10)
    cvx_exists = true;
end

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
    msg = 'nucl should be 1H, 13C or 2H\n';
    error(msg);
end
%computing the parameters from the input for scalable trajectory
theta_lower = theta(1);
theta_upper = theta(1);
RES = RES/sqrt(cos(theta_lower)^2+sin(theta_upper)^2);
theta_range = atan(sin(theta_upper)/cos(theta_lower));

% generating cone waveform
[g,k] = gencone(RES,FOV,Nint,theta_range,nucl,MAXLEN,TS,SMAX,GMAX);

% angular positions of the cones
N_cones = ceil(pi/2 * 10*FOV/RES) + 1; %multiplier 10 comes from the conversion from cm to mm
n = linspace(0,N_cones-1,N_cones-1);
polar_angles = pi*(-1/2 + n/(N_cones-1));

% scaling waveform for the stack
g_seq_x = [];
g_seq_y = [];
g_seq_z = [];
rewinderBoolean = [];
for i = 1:length(n)
    g_tmp(:,1) = cos(polar_angles(i))/cos(theta_lower)*g(:,1);
    g_tmp(:,2) = cos(polar_angles(i))/cos(theta_lower)*g(:,2);
    g_tmp(:,3) = sin(polar_angles(i))/sin(theta_upper)*g(:,3);
    
    % k-space trajectory calcualtion
    % [g] = G/cm [TS] = s [gamma] = 10kHz/T -> hence division by 100 for
    % back-conversion to MHz
    % result: [k] = 1/m
    k_traj(:,1) = cumsum(g_tmp(:,1)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
    k_traj(:,2) = cumsum(g_tmp(:,2)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
    k_traj(:,3) = cumsum(g_tmp(:,3)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio  
    
    % rotation of the cones
    phi = 4*i*pi/N_cones;
    random_angle_x = 0;
    random_angle_y = 0;

    Rot = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];    
    k_traj = k_traj*Rot;
    
    g_tmp(:,1) = [0;diff(k_traj(:,1)/gamma/TS/100)];
    g_tmp(:,2) = [0;diff(k_traj(:,2)/gamma/TS/100)];
    g_tmp(:,3) = [0;diff(k_traj(:,3)/gamma/TS/100)];
    
    % the 100 multiplier is missing on purpose and compensated in the
    % gradRewinder3D_v2
    k_tmp(1) = k_traj(end,1);
    k_tmp(2) = k_traj(end,2);
    k_tmp(3) = k_traj(end,3);
    if cvx_exists
        % rewinder calculation based on cvx solver
        [k,g_rewinder] = cvx_rewinder(k_tmp,[0 0 0],[g_tmp(end,:)]/100,[0 0 0],GMAX/100,SMAX/100,10.71e6,TS,size(g,1)/2,1e-10);
        g_rewinder = g_rewinder*100;
    else
        % manual rewinder calculation in case cvx software in not available
        [g_rewinder,k,slewRate] = gradRewinder3D([g_tmp(end,:)],[0 0 0],k_tmp,[0 0 0],GMAX,SMAX,TS,nucl);
    end
    % concatenating rawinder with the cone segment
    g_tmp_x = [g_tmp(:,1);g_rewinder(:,1)];
    g_tmp_y = [g_tmp(:,2);g_rewinder(:,2)];
    g_tmp_z = [g_tmp(:,3);g_rewinder(:,3)];
    % indexing which entries belong to the rewinder after concatenation
    rewinderBoolean_tmp = zeros(length(g_tmp_x),1);
    rewinderBoolean_tmp(end-size(g_rewinder,1)+1:end) = 1;
    rewinderBoolean = [rewinderBoolean;rewinderBoolean_tmp];
    % adding the cone segment + rewinder to the full gradient trajectory
    g_seq_x = [g_seq_x;g_tmp_x];
    g_seq_y = [g_seq_y;g_tmp_y];
    g_seq_z = [g_seq_z;g_tmp_z,i*ones(length(g_tmp_z),1)];
    % noting the length of the cone segment + rewinder for bookkeeping
    coneLengths(i) = length(g_tmp_x);
end
clear k_traj
cone_endpoints = cumsum(coneLengths);
%% plot for checking the trajectory and that the hardware constraints are met
% k-space trajectory calcualtion
% [g] = G/cm [TS] = s [gamma] = 10kHz/T -> hence division by 100 for
% back-conversipno to MHz
% result: [k] = 1/m
k_traj(:,1) = cumsum(g_seq_x*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
k_traj(:,2) = cumsum(g_seq_y*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
k_traj(:,3) = cumsum(g_seq_z(:,1)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio

% k-space trajectory
figure
plot3(k_traj(:,1),k_traj(:,2),k_traj(:,3));
title('k-space trajectory [1/m]')
xlabel('x');
ylabel('y');
zlabel('z');
fprintf('Total time = %i ms',1000*size(g_seq_x,1)*TS)

% gradient waveform
figure
T = linspace(0,1000*size(g_seq_x,1)*TS,size(g_seq_x,1));
subplot(3,1,1)
plot(T,g_seq_x);
title('Gradient timing')
xlabel('ms');
ylabel('x gradient (G/cm)');
subplot(3,1,2)
plot(T,g_seq_y);
xlabel('ms');
ylabel('y gradient (G/cm)');
subplot(3,1,3)
plot(T,g_seq_z(:,1));
xlabel('ms')
ylabel('z gradient (G/cm)');

% slew-rate
figure
subplot(3,1,1)
plot(T(1:end-1),diff(g_seq_x)/TS);
title('Slew rates (G/cm/s)')
xlabel('ms');
ylabel('x slew rate');
subplot(3,1,2)
plot(T(1:end-1),diff(g_seq_y)/TS);
xlabel('ms');
ylabel('y slew rate');
subplot(3,1,3)
plot(T(1:end-1),diff(g_seq_z(:,1))/TS);
xlabel('ms')
ylabel('z slew rate');

% velocity k-space
figure
kv_traj(:,1) = cumsum([1:length(g_seq_x)]'.*TS.*g_seq_x*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
kv_traj(:,2) = cumsum([1:length(g_seq_y)]'.*TS.*g_seq_y*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
kv_traj(:,3) = cumsum([1:length(g_seq_z)]'.*TS.*g_seq_z(:,1)*1e-2*gamma/100*1e6*TS); % converting back to 10.71MHz/T gyromagnetic-ratio
subplot(3,1,1)
plot(T,kv_traj(:,1));
title('Velocity kspace')
xlabel('ms');
ylabel('k_v position');
subplot(3,1,2)
plot(T,kv_traj(:,2));
xlabel('ms');
ylabel('k_u position');
subplot(3,1,3)
plot(T,kv_traj(:,3));
xlabel('ms')
ylabel('k_w position');

%% density calculation for the reconstruction
for i = 1:length(g_seq_x)
    g_radial = g_seq_z(i,1)/sin(polar_angles(g_seq_z(i,2)));
    g_circ = abs(g_seq_x(i) + 1i*g_seq_y(i));
    k_radial = norm(k_traj(i,:));
    density(i) = 1/norm([g_seq_x(i),g_seq_y(i),g_seq_z(i,1)])*Nint/(k_radial^2)/cos(polar_angles(g_seq_z(i,2)))*sqrt(1+(g_circ^2)/(g_radial^2));
end
density(find(isnan(density))) = max(density);
gw = ([g_seq_x,g_seq_y,g_seq_z]);






