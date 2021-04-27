function [adjMciccc,CC,PL,SW,EG] = ProbThresh_function(events,adjM,repNum,adjMmethod,downSample,lag,fs,tail,pause)

% Author RCFeord
% Updated by HSmith, Cambridge, December 2020
% Modified for validation of repNum by LBurn, Cambridge, March 2021


% This function generates a synthetic matrix of binary data
%
% Synthetic data is created using original data and shuffling events for
% each cell using 'circular' method. All synthetic data is created from
% existing events
%
% It then averages the functional connectivity established across n number
% of iterated synthetic matrices
%
% INPUTS:
%   events = binary matrix of spikes/neuronal events, columns are nodes/cells,
%            rows correspond to time points
%   adjM = adjacency matrix of real data
%   repNum = number of repetitions for the generation of synthetic datasets
%   adjMmethod = correlation method for adjacency matrix generation:
%            'tileCoef' = STTC, 'correlation','partialcorr','xcorr'
%   downSample = downsampling (if none = 0)
%   lag = lag for the correlation or STTC
%   fs = sampling frequency
%   ploton:
%       = 0 do not show figures
%       = 1 show figures
%   tail = confidence interval. Eg input 0.05 for p = 0.05 thresholding and
%          0.025 for p = 0.025 thresholding. No default set

% OUTPUTS:
%   adjMci = real adjacency matrix thresholded at "tail" confidence interval
%            of probabilistic edge weight

%% Generate synthetic data
% adjMi specifically NOT preallocated
num_frames = size(events,1);
num_nodes = size(events,2);

CC = [];
PL = [];
SW = [];
EG = [];


for n = 1:repNum
    
    % Create a matrix the same size as the real data matrix ('events')
    SynthDatBin = zeros(size(events));
    % Select points along timeseries
    locs = randi(num_frames,1,num_nodes);
    
    % Shuffle data
    for i = 1:num_nodes
        SynthDatBin(1 : end - locs(i) +1 , i) = events(locs(i) : end , i);
        SynthDatBin(end - locs(i) +2 : end , i) = events(1 : locs(i) -1 , i);
    end
    
    % Generate adjacency matrix from synthetic data
    adjMs = getAdjM_CaPop(SynthDatBin, adjMmethod, downSample, lag, fs);
    
    % Add to stack of synthetic adjacency matrices
    adjMi(:,:,n) = adjMs;
    
    
    
    if rem(n,pause)==0
        
        adjMiccc = adjMi;
        
        % Remove negatives, NaNs, and nonzero diagonals
        adjMiccc(adjMiccc<0) = 0;
        adjMiccc(isnan(adjMiccc)) = 0;
        adjMiccc = bsxfun(@times, ~eye(num_nodes), adjMiccc);
        
        % Stats test; threshold each element if >= top "tail" % of data
        adjMciccc = adjM;
        cutoff_pointccc = ceil((1 - tail) * n);
        for i = 1:size(adjMiccc,1)
            for j = 1:size(adjMiccc,2)
                Eu = sort(adjMiccc(i,j,:),'ascend');
                if Eu(cutoff_pointccc) > adjM(i,j)
                    adjMciccc(i,j) = 0;
                end
            end
        end    
        
        %PL
        L = weight_conversion(adjMciccc, 'lengths');
        D = distance_wei(L);
        PLccc = charpath(D,0,0);
        
        %CC
        CCo = clustering_coef_wu(adjMciccc);
        CCccc = mean(CCo);
        
        %SW
        SWccc = CCccc/PLccc;
        
        %Eglob
        adjM_nrm = weight_conversion(adjMciccc, 'normalize'); %rescale all edge weights to [0,1]
        EGccc = efficiency_wei(adjM_nrm); %normalization = binarization + ??
        
        CC = [CC CCccc];
        PL = [PL PLccc];
        SW = [SW SWccc];
        EG = [EG EGccc];
        
    end
    
end
end
