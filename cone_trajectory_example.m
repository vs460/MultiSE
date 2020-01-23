% The script generates the gradient waveforms for the cone trajectory and 
% the segmented readout presented in the paper:
% 
%   Somai V, Wright AJ, Fala M, Hesse F, Brindle KM. A multi spin echo pulse 
%   sequence with optimized excitation pulses and a 3D cone readout for 
%   hyperpolarized 13C imaging. Magn. Reson. Med 2020. 
% 
%
% Vencel Somai 2020. -> vs460@cam.ac.uk


clear all
close all

nucl          = '13C';           % nucleus of interest
FOV           = 3.2;             % field of view [cm], clinical example: 16
RES           = 2;               % resolution [mm],    clinical example: 5
imageDim      = 10*FOV/RES;      % image dimension = FOV/RES
Nint          = 1;               % number of interleaves

%% generating k-space trajectory
[g,k,dens,N_cones,coneLengths,rewinderBoolean] = ...                           
    coneTrajDesign(RES, FOV/2,Nint, [0+1e-5 pi/2-1e-5],nucl);
cone_endpoints = cumsum(coneLengths);
cone_endpoints = [0,cone_endpoints];

%% generating the segmented readout and write the waveforms into file (VnmrJ format)
path = uigetdir();               % specify folder to save the waveform files to
for cc = 1:floor((N_cones-1)/2)
    % X grad
    Xgradamp = 0;
    fileID = fopen(fullfile(path,strcat('vs_cones_',nucl,'_',num2str(imageDim),'X_dual_clin',num2str(cc),'.GRD')),'w');
    for n = cone_endpoints(cc+1):-1:cone_endpoints(cc)+1
        A = -g(n,1)*32767/max(abs(g(:,1)));
        if abs(g(n,1)) > Xgradamp
            Xgradamp = abs(g(n,1));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    for n = cone_endpoints(N_cones-cc)+1:cone_endpoints(N_cones-cc+1)
        A = g(n,1)*32767/max(abs(g(:,1)));
        if abs(g(n,1)) > Xgradamp
            Xgradamp = abs(g(n,1));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    fclose(fileID);
    
    % Y grad
    Ygradamp = 0;
    fileID = fopen(fullfile(path,strcat('vs_cones_',nucl,'_',num2str(imageDim),'Y_dual_clin',num2str(cc),'.GRD')),'w');
    for n = cone_endpoints(cc+1):-1:cone_endpoints(cc)+1
        A = -g(n,2)*32767/max(abs(g(:,2)));
        if abs(g(n,2)) > Ygradamp
            Ygradamp = abs(g(n,2));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    for n = cone_endpoints(N_cones-cc)+1:cone_endpoints(N_cones-cc+1)
        A = g(n,2)*32767/max(abs(g(:,2)));
        if abs(g(n,2)) > Ygradamp
            Ygradamp = abs(g(n,2));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    fclose(fileID);
    
    % Z grad
    Zgradamp = 0;
    fileID = fopen(fullfile(path,strcat('vs_cones_',nucl,'_',num2str(imageDim),'Z_dual_clin',num2str(cc),'.GRD')),'w');
    for n = cone_endpoints(cc+1):-1:cone_endpoints(cc)+1
        A = -g(n,3)*32767/max(abs(g(:,3)));
        if abs(g(n,3)) > Zgradamp
            Zgradamp = abs(g(n,3));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    for n = cone_endpoints(N_cones-cc)+1:cone_endpoints(N_cones-cc+1)
        A = g(n,3)*32767/max(abs(g(:,3)));
        if abs(g(n,3)) > Zgradamp
            Zgradamp = abs(g(n,3));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    fclose(fileID);
end
if mod(N_cones-1,2)
    cc = cc + 1;
    
    % X grad
    Xgradamp = 0;
    fileID = fopen(fullfile(path,strcat('vs_cones_',nucl,'_',num2str(imageDim),'X_dual_clin',num2str(cc),'.GRD')),'w');
    for n = cone_endpoints(cc)+1:cone_endpoints(cc+1)
        A = g(n,1)*32767/max(abs(g(:,1)));
        if abs(g(n,1)) > Xgradamp
            Xgradamp = abs(g(n,1));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    fclose(fileID);
    
    % Y grad
    Ygradamp = 0;
    fileID = fopen(fullfile(path,strcat('vs_cones_',nucl,'_',num2str(imageDim),'Y_dual_clin',num2str(cc),'.GRD')),'w');
    for n = cone_endpoints(cc)+1:cone_endpoints(cc+1)
        A = g(n,2)*32767/max(abs(g(:,2)));
        if abs(g(n,2)) > Ygradamp
            Ygradamp = abs(g(n,2));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    fclose(fileID);
    
    % Z grad
    Zgradamp = 0;
    fileID = fopen(fullfile(path,strcat('vs_cones_',nucl,'_',num2str(imageDim),'Z_dual_clin',num2str(cc),'.GRD')),'w');
    for n = cone_endpoints(cc)+1:cone_endpoints(cc+1)
        A = g(n,3)*32767/max(abs(g(:,3)));
        if abs(g(n,3)) > Zgradamp
            Zgradamp = abs(g(n,3));
        end
        fprintf(fileID,'%6.0f 1.0\n',A);
    end
    fclose(fileID);
end

