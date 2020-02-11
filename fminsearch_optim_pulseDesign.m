% The script generates the excitation RF waveforms for the multi spin echo 
% pulse sequence presented in the paper:
% 
%   Somai V, Wright AJ, Fala M, Hesse F, Brindle KM. A multi spin echo pulse 
%   sequence with optimized excitation pulses and a 3D cone readout for 
%   hyperpolarized 13C imaging. Magn. Reson. Med 2020. 
%
% Vencel Somai 2020. -> vs460@cam.ac.uk

clear all
close all

%% set parameters
N     = 402;                                % number of initial filter points
time  = 11.256e-3                           % total duration of the pulse [s]
f     = 1/time/*N;                      % physical smapling bandwidth
om_s  = 400;                                % stop-band endge [Hz]
om_p  = 200;                                % pass-band endge [Hz]
F     = [0 om_p om_s f]/f;                  % array of relative frequencies for the fir design funtion
Amp   = [1/sqrt(2) 1/sqrt(2) 0 0];          % array of amplitudes of each freqeuncy region
W     = [1 1];                              % weight of each frequency region
BW    = (om_p + om_s)/2/pi*f;               % normalized pass-band width
gamma = 10.71e6;                            % 13C gyromagnetic ratio [H/T]
TS    = 1/f;                                % sampling time of the pulse [s]
T1    = 50;                                 % T1 relaxation time for [1-13C]lac
T2    = 0.2;                                % T2 relaxation time for [1-13C]lac

%% polynomial generation  and iSLR
beta = firls(N,F,Amp,W);                    % generating linear phase filter
alpha = gen_alpha(beta);                    % generating matching minimum-phase alpha 
rf = iSLR(alpha,beta,gamma,TS);             % inverse-SLR transform
rf_init = -imag(rf);                        
tau = time/N;                               % rounding it to the nearest multiple of 4e-6s
time = length(rf_init)*ceil(tau/4e-6)*4e-6; % 
rf_init = [rf_init;...                      % reshaping the complex pulse to an array of the real and imag part for the optimization
           zeros(1,length(rf_init))];
%% refine it with optimization
Nf = 700;                                   % number of frequency points where the result is evaluated in the cost fcn
ideal_profile = zeros(1,Nf);               
ideal_profile(501:end) = sin(45/180*pi);    % design goal specification
options = optimset('MaxFunEvals',1e8,'MaxIter',1e8,'TolFun',1e-6,'TolX',1e-6);
pulse = fminsearch(@(pulse) fmin_pulse_design_cost(pulse,ideal_profile,Nf,time,N,T1,T2),rf_init,options);
pulse = pulse(1,:) + 1i*pulse(2,:);         % concatenating the real and imag part to a complex waveform


%% simulate result
df = [linspace(-2000,2000,800)];            % frequency range where the result is evaluated
[mx,my,mz] = bloch1(pulse,0,time/N,15,4,df,0,0);
% transverse component
figure
plot(df,abs(mx+1i*my));
title('Bloch simulated m_{xy} profile')
% longitudinal component
figure
plot(df,mz)
title('Bloch simulated m_z profile')



