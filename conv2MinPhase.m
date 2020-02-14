function [beta] = conv2MinPhase(beta)

% The function converts the linear phase FIR filter to maximal phase for
% minumal phase RF pulse by flipping the roots located outside the unit
% circle.

N            = length(beta);                    % original filter length  
b_fir        = interp1(beta,linspace(1,N,200)); % downsample to speed up the calculations
%% symbolic evaluation of the roots and zero-plot
digits(16);                                     % desired precision for the symbolic calculations
b_fir_symb   = poly2sym(b_fir);                 % conversion to symbolic polynomial for increased precision
b_zeros_symb = vpasolve(b_fir_symb);            % solving for the roots
b_zeros      = double(b_zeros_symb);            % converting back from symbolic to double
zplane(b_zeros)                                 % zero-pole plot
title('Zero-pole plot of the original filter')
pause()
close all
%% flip the pass-band zeros located outside the unit circle
zero_indeces               = find((abs(b_zeros)-1)>0.01);         % searching for roots outside the unit circe
b_zeros_symb(zero_indeces) = 1./conj(b_zeros_symb(zero_indeces)); % zero flip
zplane(double(b_zeros_symb));                                     % zero-pole plot
title('Zero-pole plot after flipping the roots')
pause()
close all
%% create the final polynomial from the zeros and zero-plot to 
v         = [b_zeros_symb];                                       
x         = sym('x');
yourpoly  = expand(prod(x-v));                                    % building the final polynomial from the roots
b_fir_tmp = sym2poly(yourpoly);                                   % convertig back from symbolic 
b_fir     = b_fir_tmp;
b_zeros   = roots(b_fir);                                         % check the result once again
zplane(b_fir);
title('Zero-pole plot of the result to check numerical stability of the root search')
pause()
close all

%% interpolation back to the original length
beta      = interp1(b_fir,linspace(1,length(b_fir),N));
end

