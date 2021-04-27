% Script to validate network null model stability and quality.
% 
% 
% REQUIREMENTS: latmio_valid_function.m;
%               randmio_valid_function.m;
%               AdjMProbabilisticThresh_v3.m; 
%               getAdjM_CaPop.m
% 
% 
% LBurn, November 2020


%% start-up
HomeDir = '/Users/lance/Documents/MATLAB/My stuff/Null Models/InputData';
cd(HomeDir);

load('170510c4fr2CaPopProj.mat');
% sampling frequency (for calcium imaging 2.7-3 Hz -> approx 3 frames per second)
fs = str2double(Info.FrRate)/1000;

%% generate the adjacency matrix from the binary data

% adjacency matrix parameters:
method = 'tileCoef'; % 'tileCoef' = STTC, 'correlation','partialcorr','xcorr'
downSample = 0; % this does not apply to current dataset as temporal resolution 
% too low. Keep at 0 for now.
lag = 1; % time window for correlation in seconds
ploton = 0;
tail = 0.05;
repNum = 1500;

adjM = getAdjM_CaPop(BinEvents, method, downSample, fs, lag);

% remove negatives
adjM(adjM<0) = 0;
adjM(isnan(adjM)) = 0;
for i = 1:length(adjM)
    adjM(i,i) = 0;
end

% threshold AdjM (probabilistic)
[W_thr] = AdjMProbabilisticThresh_v3(BinEvents,adjM,repNum,method,downSample,lag,fs,ploton,tail);


%% Obtain data for validation of iterations

%   **** Inputs ****
data = 10;              % How many times you want to run through the script (to get a mean and SD; e.g., 10).
ITERATIONS = 10000;     % How many iterations to include per initiation (e.g., 10,000).
nullmodel = 'latt';     % What null model do you want to validate (rand, latt).  
pause = 10;             % How often do you want to pull out the metrics (e.g., every 10 iterations).

t=1;
CC_all=[];
PL_all=[];
SW_all=[];
EG_all=[];

Z = pdist(W_thr);
D = squareform(Z);

if strcmp(nullmodel,'rand')
    for ii=1:data
        [R,CC,PL,SW,EG]=randmio_valid_function(W_thr, ITERATIONS, pause);
        CC_all(t,:) = CC;
        PL_all(t,:) = PL;
        SW_all(t,:) = SW;
        EG_all(t,:) = EG;
        t=t+1;
    end
end

if strcmp(nullmodel,'latt')
    for ii=1:data
        [R,CC,PL,SW,EG]=latmio_valid_function(W_thr, ITERATIONS, D, pause);
        CC_all(t,:) = CC;
        PL_all(t,:) = PL;
        SW_all(t,:) = SW;
        EG_all(t,:) = EG;
        t=t+1;
    end
end


%% Real metrics

% real CC
K=sum(W_thr~=0,2);
cyc3=diag((W_thr.^(1/3))^3);
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
Cc=cyc3./(K.*(K-1));        %real clustering coefficient
CCreal = mean(Cc);

%PL
Ln = weight_conversion(W_thr, 'lengths');
D = distance_wei(Ln);
PLreal = charpath(D,0,0);

%SW
SWreal = CCreal/PLreal;

%EG
adjM_nrm = weight_conversion(W_thr, 'normalize'); %rescale all edge weights to [0,1]
Eglobreal = efficiency_wei(adjM_nrm); %normalization = binarization + ??


%% Plotting



met_mean = mean(CC_all);
met_SD = std(CC_all);

x = 1:10:ITERATIONS;

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

yline(CCreal, '-k', 'Empirical network value');
%xlim([? ?]);
%ylim([0.03 0.13]);
xline(3000, '--r', 'Chosen ITER value');
xlabel('Iteration', 'FontSize',15);
ylabel('Clustering coefficient', 'FontSize',15);
aesthetics();


%% Save
cd('/Users/lance/Documents/MATLAB/My stuff/Null Models/OutputData/Iteration Validation');

save('Variables_latt');
saveas(F1,'SW_rand_line.fig');