function [rf] = iSLR(alpha,beta,gamma,TS)
% Performs inverse SLR transform based on the paper:
%
%   Pauly J, Le Roux P, Nishimura D, Macovski A. Parameter relations for the
%   Shinnar-Le Roux selective excitation pulse design algorithm [NMR imaging].
%   IEEE Trans Med Imaging. 1991;10(1):53-65
% 
%   Input:
%       alpha = coefficients of the alpha polynomial  
%       beta  = coefficients of the beta polynomial  
%       gamma = gyromagnetic ratio [Hz/T]
%       TS    = RF pulse sampling time [s]
%   Output:
%       rf    = RF waveform
% 

for i = length(alpha):-1:1 
    %% Calculate the ith element of the RF waveform
    phi = 2*atan2(abs(beta(1)),abs(alpha(1)));
    theta = angle(-1i*beta(1)/alpha(1));
    rf(i) = 1/gamma/TS*phi*exp(1i*theta);
    
    %% Do the inverse recursion 
    C = cos(abs(rf(i)/2));
    S = 1i*exp(1i*angle(rf(i)))*sin(abs(rf(i)/2));
    Nalpha = length(alpha);
    alpha_tmp = C*alpha + conj(S)*beta;
    beta_tmp = conv([-S zeros(1,Nalpha)],alpha) + conv([C zeros(1,Nalpha)],beta);
    alpha = alpha_tmp(1:Nalpha-1);
    beta = beta_tmp(2:Nalpha);
end

end






