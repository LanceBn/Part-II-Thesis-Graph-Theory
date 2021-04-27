% Script to validate the repNum needed for probabilistic thresholding.

% INPUTS:
%
%   Data    = How many times you repeat the validation (to obtain S.D.)
%   repNum  = Repetition number, the maximum value you want to validate to
%   pause   = Measure network metrics every 'n' repNum
%   tail    = Set the p-value for probabilistic thresholding

% OUTPUT
%   A set of network metrics (topologies) recorded from the post-threshold
%   matrix as it is progressively thresholded with increasing repNum. Plot
%   these metrics against repNum to visualise how the process stabilises
%   with increasing repNum.

% REQUIREMENTS
%   ProbThresh_valudation.m function is required, as it is called to within
%   this script.
%   getAdjM_CaPop.m is needed to create the adjacency matrix.


% Lance Burn, March 2021

%% start-up
HomeDir = '/Users/lance/Documents/MATLAB/My stuff/Probabilistic Thresholding';
cd(HomeDir);
cd('InputData');

load('170510c4fr2CaPopProj.mat');
cd(HomeDir);
% sampling frequency (for calcium imaging 2.7-3 Hz -> approx 3 frames per second)
fs = str2double(Info.FrRate)/1000;

%% generate the adjacency matrix from the binary data

% adjacency matrix parameters:
method = 'tileCoef'; % 'tileCoef' = STTC, 'correlation','partialcorr','xcorr'
downSample = 0; % this does not apply to current dataset as temporal resolution 
% too low. Keep at 0 for now.
lag = 1; % time window for correlation in seconds

adjM = getAdjM_CaPop(BinEvents, method, downSample, fs, lag);

% remove negatives
adjM(adjM<0) = 0;
adjM(isnan(adjM)) = 0;
for i = 1:length(adjM)
    adjM(i,i) = 0;
end



%% Obtain data for validation of repNum

%   **** Inputs ****
data = 10;              % how many times you want to create a null model of 'ITERATIONS' number of iters/
repNum = 10000;     % How many iterations to include per run (note - we sample every 10 iterations).
pause = 10;
tail = 0.05;

t=1;
CC_all=[];
PL_all=[];
SW_all=[];
EG_all=[];


% % Real CC (n = 0)
% CCo = clustering_coef_wu(W_thr);
% CCom = mean(CCo);
% 
% %Real PL
% ...
%     
% %Real SW
% ...
    

for ii=1:data
    [adjMci,CC,PL,SW,EG] = ProbThresh_function(BinEvents,adjM,repNum,method,downSample,lag,fs,tail,pause);
    CC_all(t,:) = CC;
    PL_all(t,:) = PL;
    SW_all(t,:) = SW;
    EG_all(t,:) = EG;
    t=t+1;
end



%% Plotting

met_mean = mean(CC_all);
met_SD = std(CC_all);

x = 1:pause:repNum;

% plot
F1 = figure;
SD_curveU = met_mean + met_SD;
SD_curveL = met_mean - met_SD;
Xf =[x,fliplr(x)];                                      % create continuous x value array for plotting
Yf =[SD_curveU,fliplr(SD_curveL)];                      % create y values for out and then back
h = fill(Xf,Yf,[0, 0.4470, 0.7410],'edgecolor','none'); 
set(h,'facealpha',0.3)                                  % 0 (invisible) and 1 (opaque). 
hold on
plot(x,met_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',3)

% yline(wo, '-b', 'Actual network small-worldness, w');
%xlim([? ?]);
%ylim([0.03 0.13]);
xlabel('repNum', 'FontSize',15);
ylabel('Clustering coefficient', 'FontSize',13);
% title('');
aesthetics();


%% Save
cd('/Users/lance/Documents/MATLAB/My stuff/Null Models/OutputData/Iteration Validation');

save('repNum10000x10');
saveas(F1,'CC.fig');
