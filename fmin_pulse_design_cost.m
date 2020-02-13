function [ result ] = fmin_pulse_design_cost(pulse,goal,Nf,time,Np,T1,T2,df,inhomHz,target_metab)

% Cost function for the optimized excitation pulse design used in paper:  
% 
%   Somai V, Wright AJ, Fala M, Hesse F, Brindle KM. A multi spin echo pulse 
%   sequence with optimized excitation pulses and a 3D cone readout for 
%   hyperpolarized 13C imaging. Magn. Reson. Med 2020.  
% 
% ****************************************************************************
% !!!Not recommended for routine pulse design due to computatinal cost issues,
% only as a targeted tool!!!
% ****************************************************************************
%
%   Input:
%       pulse = 2xN array with the real and the imaginary part of the RF
%       waveform
%       goal  = the ideal profile 
%       Nf    = number of frequency points where the result is evaluated
%       time  = duration of the pulse [s]
%       Np    = number of (complex) pulse points
%       T1    = T1 relaxation time for Bloch-simulation [s]
%       T2    = T2 relaxation time for Bloch-simulation [s]
%   Output:
%       result = cost generated by a given solution (RF waveform) 
% 
% Vencel Somai 2020. -> vs460@cam.ac.uk 

offset      = strcmp(target_metab,'pyr')*916;   % offset frequency between lac and pyr     
weights     = ones(size(df));                   % weights indexing where the frequency profile is of interest
weights_lac = ones(size(df));                   % weights indexing the lactate frequency band
weights_pyr = ones(size(df));                   % weights indexing the pyruvate frequency band
weights_ala = ones(size(df));                   % weights indexing the lactate alanine band
weights_PH  = ones(size(df));                   % weights indexing the pyruvate-hydrate frequency band
weights_lac(df<offset-inhomHz)       = 0;       % specifying the lactate band, the weights are 0 otherwise
weights_lac(df>offset+inhomHz)       = 0;       % specifying the lactate band, the weights are 0 otherwise
weights_pyr(df<offset-916-3*inhomHz) = 0;       % specifying the pyruvate band, the weights are 0 otherwise
weights_pyr(df>offset-916+3*inhomHz) = 0;       % specifying the pyruvate band, the weights are 0 otherwise
weights_ala(df<offset-487-inhomHz)   = 0;       % specifying the alanine band, the weights are 0 otherwise
weights_ala(df>offset-487+inhomHz)   = 0;       % specifying the alanine band, the weights are 0 otherwise
weights_PH(df<offset-300-inhomHz)    = 0;       % specifying the pyruvate-hydrate, the weights are 0 otherwise
weights_PH(df>offset-300+inhomHz)    = 0;       % specifying the pyruvate-hydrate band, the weights are 0 otherwise
weights = weights_lac + weights_pyr + weights_ala + weights_PH;

%% Projection to a set of RF pulses with given smoothness
% Not proved mathemtically but helps to speed up convergence to a solution
% and the measured frequency profile mathes the Bloch-simulations much 
% better if the RF waveform is somewhat smooth.
pulse_tmp = smooth(pulse(1,:),9)' + 1i*smooth(pulse(2,:),9)';

% evaluate the result at a range of B1 amplitudes
% for the pre-clinical volume coil used 0.5 - 2 range is reasonable
B1 = 0.5:0.2:2;
    pulse = B1'*pulse_tmp;
    mx = zeros(length(B1),length(df));my = zeros(length(B1),length(df));mz = zeros(length(B1),length(df));
    parfor i = 1:length(B1)
        [mx(i,:),my(i,:),mz(i,:)] = bloch1(pulse(i,:),0,time/Np,T1,T2,df,0,[0 0]);
        while sum(isnan(mx(i,:))+isnan(my(i,:))+isnan(mz(i,:)))>0
            [mx(i,:),my(i,:),mz(i,:)] = bloch1(pulse(i,:),0,time/Np,T1,T2,df,0,[0 0]);
        end
    end
    if strcmp(target_metab,'lac')
        result = sum(sum(repmat(weights,[length(B1),1]).*(abs(mx+1i*my) - repmat(goal,[length(B1),1])).^2 + ...
                 repmat(2*weights_lac,[length(B1),1]).*(angle(mx+1i*my) - angle(repmat(goal,[length(B1),1]))).^2))+...
                 1e6*(max(max(abs(pulse)))>maxB1)
    else
        result = sum(sum(repmat(weights,[length(B1),1]).*(abs(mx+1i*my) - repmat(goal,[length(B1),1])).^2 + ...
                 repmat(2*weights_pyr,[length(B1),1]).*(angle(mx+1i*my) - angle(repmat(goal,[length(B1),1]))).^2))+...
                 1e6*(max(max(abs(pulse)))>maxB1)
    end        
end

