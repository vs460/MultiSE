function [alpha] = gen_alpha(beta)
% Calculates the alpha polynomial from given beta polynomial

N = length(beta);
beta_zf = [beta zeros(1,2^ceil(log2(N)) - N)]; % zero-fill b to the next power of 2 for fft
B = fft(beta_zf);
A_mag = sqrt(1 - B.*conj(B));                  % apply magnitude constraint 
A = exp(ifft(hilbert(fft(log(A_mag)))));            % apply formula for minimum phase alpha
alpha = fft(A)/(2^ceil(log2(N)));              % normalization
alpha = real(alpha(1:N));                      % use only the original points

end

