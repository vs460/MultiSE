function [k,g] = cvx_rewinder(k0,kend,g0,gend,GMAX,SMAX,gamma,TS,N_upb,eps)

% The function attempts to find time optimal rewinder design based on the solver
% "CVX: Matlab Software for Disciplined Convex Programming":
%   http://cvxr.com/cvx/
%   email: info@cvxr.com
% 
%   [1] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined 
%   convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.
%
%   [2] Michael Grant and Stephen Boyd. Graph implementations for nonsmooth 
%   convex programs, Recent Advances in Learning and Control 
%   (a tribute to M. Vidyasagar), V. Blondel, S. Boyd, and H. Kimura, editors,
%   pages 95-110, Lecture Notes in Control and Information Sciences, Springer, 
%   2008. http://stanford.edu/~boyd/graph_dcp.html.
%
% Input:
%
%   k0    = initial point in k-space [1/m]
%   kend  = desired endpoint in k-space [1/m]
%   g0    = initial gradient values [T/m]
%   gend  = desired gradient values at the end of the waveform [T/m]
%   GMAX  = maximal gradient amplitude [T/m]
%   SMAX  = maximal slew-rate [T/m/s]
%   gamma = gyromagnetic ratio (\gamma/(2\pi)) [Hz/T]
%   TS    = sampling time of the gradient waveform [s]
%   N_upb = upper bound of the length of the wavefom for the bisectioning
%   eps   = tolerance with respect to the resulted k-space displacement
% 
% Output:
%
%   k     = resulted k-space trajectory [1/m]
%   g     = resulted gradient waveform [T/m]
%
% **************************************************************************
% !!!gradient axes are treated separately with respect to hardware limits!!!
% The paper:
%
%
%
% Take care if used for general purpose, not tested extensively! 
% No responsibilty taken!
% **************************************************************************
%
% Vencel Somai 2020. -> vs460@cam.ac.uk


%% 
addpath('C:\Users\Somai01\Documents\PulseDesign\cvx')
N_lwb = 0;          % lower bound of the length of the waveform for bisectioning
N     = round(N_upb/2);    % initial value for the length of the waveform

%% 
repetitive_flag = 0;
exit_flag = false;
while (~(exit_flag))
    D   = toeplitz([-1/TS 1/TS zeros(1,N+1)]',zeros(1,N+3))'; % numerical differentiation mtx
    % ********************************* start solver *************************** %
    cvx_begin
    cvx_precision low
    variable g(N,3);
    dk1 = sum(g(:,1))*gamma*TS;             % resulted displacement along kx axis
    dk2 = sum(g(:,2))*gamma*TS;             % resulted displacement along ky axis
    dk3 = sum(g(:,3))*gamma*TS;             % resulted displacement along kz axis
    s   = D*[g0;g;gend;gend];
    s = s(1:end-1);
    minimize norm(s(:,1),2) + norm(s(:,2),2) + norm(s(:,3),2) % to make use of full slew-rate
    subject to % ********************** constraints **************************** %
    s               <= SMAX;                % upper bound for the slew-rate
    s               >= -SMAX;               % lower bound for the slew-rate
    g               <= GMAX;                % upper bound for the gradient amplitude
    g               >= -GMAX;               % lower bound for the gradient amplitude
    [dk1,dk2,dk3]   <= [kend-k0]*(1+eps);   % displacement in the k=spce should match the desired within tolerance
    [dk1,dk2,dk3]   >= [kend-k0]*(1-eps);   % displacement in the k=spce should match the desired within tolerance
    cvx_end
    % **********************************end solver ***************************** %
    
    % check whether the resulted waveformn is feasible
    if (isnan(cvx_optval) || isinf(cvx_optval) || max(max(abs(s)))>SMAX*(1+eps) || max(max(abs(g)))>GMAX*(1+eps) || sum(sum(isinf(g)))>0 || sum(sum(isnan(g)))>0 || sum(sum(isinf(s)))>0 || sum(sum(isnan(s)))>0)
        N_lwb_old = N_lwb;
        N_lwb = N;
        N = ceil((N_upb - N_lwb)/2 + N_lwb);
    else
        N_old = N;N_upb_old = N_upb;N_lwb_old = N_lwb;
        N_upb = N;
        N = floor((N_upb - N_lwb)/2 + N_lwb);
        if (N_lwb == N_lwb_old && N_upb == N_upb_old) % if there is no more change in the interval exit the loop
            exit_flag = true;
        end
    end
end

% calculate the resulted k-sapce segment
k = [cumsum(g(:,1))*gamma*TS,cumsum(g(:,2))*gamma*TS,cumsum(g(:,3))*gamma*TS];
%%  plot
T = linspace(0,size(g,1)*TS,size(g,1))*1e3;
% *********************** resulted k-space trajectory *************************
figure
plot3(k(:,1),k(:,2),k(:,3))
sgtitle('k-Space trajectory [1/m]')
% *****************************************************************************

% *********************** resulted gradient waveform **************************
figure
sgtitle('Gradient waveforms')
subplot(3,1,1)
plot(T,g(:,1))
xlabel('Time [ms]')
ylabel('Gx [G/cm]')
subplot(3,1,2)
plot(T,g(:,2))
xlabel('Time [ms]')
ylabel('Gy [G/cm]')
subplot(3,1,3)
plot(T,g(:,3))
xlabel('Time [ms]')
ylabel('Gz [G/cm]')
% *****************************************************************************

% ************************** resulted slew-rate *******************************
figure
sgtitle('Slew-rates')
subplot(3,1,1)
plot(T(1:end-1),diff(g(:,1))/TS)
xlabel('Time [ms]')
ylabel('Sx [G/cm/ms]')
subplot(3,1,2)
plot(T(1:end-1),diff(g(:,2))/TS)
xlabel('Time [ms]')
ylabel('Sy [G/cm/ms]')
subplot(3,1,3)
plot(T(1:end-1),diff(g(:,3))/TS)
xlabel('Time [ms]')
ylabel('Sz [G/cm/ms]')
% *****************************************************************************

end
