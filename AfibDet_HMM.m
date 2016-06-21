% PROJECT:  AfibDet 
%------------------
% Title:    Analysis of features of qtdb
% Autor:    Barbara Jesacher
% Date:     6. Juni 2016
% Version:  1.0
%----------------------------------------
% read signals from QT database
% HMM
% plot 
% functions used:
%   - readFeaturesFromDATA.m


%% LOAD DATA
clear
close all

addpath(genpath('/Users/barbara/Documents/MATLAB/HMMall'));
addpath(genpath('/Users/barbara/Dropbox/BFH/MatlabFunctions'))
addpath(genpath('/Users/barbara/Dropbox/BFH/MatlabScripts/'))

path = '/Users/barbara/Dropbox/BFH/HMM_QTdb/';

load(strcat(path, 'TrainData/TrainData.mat'))

%% set plot parameters
line1 = {'.b','.g','.r','.c','.k','.m'};
% line2 = {'--o','--s','--d','-->','--<','--v'};
% marker1 = {'o','s','d','>','<','v','^','h'};
% colorLabel = {'go','rs','md','c+','c*'};
color = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[0 0 0.65]};


%% Initialise 
% Q = 6;                          %hidden states
% QL = [11 20 12 315 17 60];            %mean length of each state
% Q = 3;
% QL = [35 325 65];
Q = 4;
QL = [10 25 320 65];            %mean length of each state
% 
selFeat = [1, 2, 3, 4, 5, 6];          %not more features than states!!
Fs = 250;
tol = 1e-9;
nbMaxStep = 1000;
%nbLearn = 30;
% nbWaves = 10; 
%ql = round(QL*Fs/1000);

% waveTrain.feat = readFeaturesFromDAtA(Pwave, QRScomp, Twave, idle, nbWaves, 1);
% 
% for jj = 1 : Q
%     histData = [];
%     ind1 = cumsum([0, QL])./sum(QL);
%     
%     for kk = 1 : nbWaves
%         sl = length(waveTrain.feat{1,kk});
%         ind = round(ind1(jj)*sl) + 1:round(ind1(jj+1)*sl);
%         histData = [histData; waveTrain.feat{1, kk}(selFeat, ind)'];
%     end
% 
%     HMM.meanInit(:,jj) = mean(histData)';
%     HMM.sigmaInit(:,:,jj) = cov(histData);
%end
%     HMM.mixInit = ones(Q, 1);
%     HMM.priorInit = [1; zeros(Q-1,1)];
%     if Q > 4
%         HMM.trInit = diag(QL./(QL+1)) + circshift(diag(1./(QL+1)), [0,1]);      
%     elseif Q == 4
%         HMM.trInit = diag(QL./sum(QL));%diag(QL./(QL+1));
%         HMM.trInit(1,2:4) = 1./(QL(1)+1) / 3;
%         HMM.trInit(2, 1) = 1./(QL(2)+1);
%         HMM.trInit(3, 1) = 1./(QL(3)+1);
%         HMM.trInit(4, 1) = 1./(QL(4)+1);
%     elseif Q == 3
%         HMM.trInit = diag(QL./(QL+1)) + circshift(diag(1./(QL+1)), [0,1]); 
%         
%     end
% 
% histData = [];

%%
nbBeats = 60;
nSig = 1; % [1, 3];

trainSequence = readFeaturesFromDAtA(Pwave, QRScomp, Twave, idle, nbBeats, nSig);
reshapedData = [];
for kk = 1 : nbBeats
    waveTrain.seq{kk} = trainSequence{1, kk}(selFeat, :);
    reshapedData = [reshapedData, trainSequence{1, kk}(selFeat, :)];
end
trainSequence = [];

%% Fit a mixture of M=2 Gaussians for each of the states using K-means:
M = 6;

[mu0, Sigma0] = mixgauss_init(Q*M, reshapedData, 'full');
HMM.meanInit = reshape(mu0, [length(selFeat) Q M]);
HMM.sigmaInit = reshape(Sigma0, [length(selFeat) length(selFeat) Q M]);
HMM.mixInit = ones(Q,M);
HMM.priorInit = [1; zeros(Q-1,1)];
if Q > 4
    HMM.trInit = diag(QL./(QL+1)) + circshift(diag(1./(QL+1)), [0,1]);      
elseif Q == 4
    HMM.trInit = diag(QL./sum(QL));%diag(QL./(QL+1));
    HMM.trInit(1,2:4) = 1./(QL(1)+1) / 3;
    HMM.trInit(2, 1) = 1./(QL(2)+1);
    HMM.trInit(3, 1) = 1./(QL(3)+1);
    HMM.trInit(4, 1) = 1./(QL(4)+1);
elseif Q == 3
    HMM.trInit = diag(QL./(QL+1)) + circshift(diag(1./(QL+1)), [0,1]); 
end

%% Train 
[~, HMM.priorProb, HMM.trProb, HMM.mean, HMM.sigma, HMM.mix] = mhmm_em(waveTrain.seq, ...
                HMM.priorInit, HMM.trInit, HMM.meanInit, HMM.sigmaInit,...
                HMM.mixInit, 'thresh', tol, 'max_iter', nbMaxStep, 'verbose',...
                true, 'cov_type', 'full');

 
                    
%% HMM Visualize
for jj = 1 : length(waveTrain.seq)
    seq1 = mat2cell(waveTrain.seq{1, jj}, [size(waveTrain.seq{1, jj}, 1), []]);
    obslik = mixgauss_prob(seq1{1}, HMM.mean, HMM.sigma, HMM.mix);
    waveTrain.seqViterbi{jj} = viterbi_path(HMM.priorProb, HMM.trProb, obslik)';
end
obslik = [];


%% Plot 
feat1 = 5;
feat2 = 6;
f1 = find(selFeat == feat1);
f2 = find(selFeat == feat2);

points = [];
states = [];
for jj = 1 : nbBeats
    states1 = waveTrain.seqViterbi{1, jj};
    points1 = [waveTrain.seq{1, jj}(f1, :)', waveTrain.seq{1, jj}(f2, :)'];
    points = [points; points1];
    states = [states; states1];
end

figure
subplot(121)
hold on
plotFeatures(idle, Pwave, QRScomp,Twave, feat1, feat2, Q, nbBeats, line1)
hold off
subplot(122)
hold on
for ii = 1 : length(points)
   plot(points(ii, 1), points(ii, 2), line1{states(ii)})
end
hold off
title('Features Viterbi - left: true path - right: most probable path')
xlabel('du_dt')
ylabel('count_dt')
if Q == 4
    legend('idle', 'P-wave', 'QRS-complex', 'T-wave')
    str = {'QL = [30 20 315 60], selFeat = [1, 2, 3, 4, 5, 6], nbMaxStep = 30, nbWaves = 50'};
elseif Q == 6
    legend('idle1', 'P-wave', 'idle2', 'QRS-complex', 'idle3', 'T-wave')
    str = {'QL == [11 20 12 315 17 60], selFeat = [3,4,5], nbMaxStep = 30, nbWaves = 50'};
elseif Q == 3
    legend('idle1', 'P-wave', 'idle2', 'QRS-complex', 'idle3', 'T-wave')
    str = {'QL == [30 330 70], selFeat = [1,2,3,4,5,6], nbMaxStep = 30, nbWaves = 50'};
end
text(20, 2.4, str)



%% Plot
load(strcat(path, 'sel103_Data.mat'))
s = Data{1,2}(1).dataAsynch.AsynchSignals.signalAsynch;
t = Data{1,2}(1).dataAsynch.AsynchSignals.timeAsynch;

figure;
title('Segmented ECG signal')
plotECGSegmented(Q, t, s, waveTrain.seqViterbi, line1)
if Q == 4
    str = {'QL == [30 20 315 60]'};
elseif Q == 6
    str = {'QL == [11 20 12 315 17 60]'};
elseif Q == 3
    str = {'QL == [30 330 70]'};
end
text(20, 2.4, str)


%% TEST
nBeats = [121, 180];
nbSig = 1;

testSequence = readFeaturesFromDAtA(Pwave, QRScomp, Twave, idle, nBeats, nbSig);

for kk = 1 : length(testSequence)
    waveTest.seq{kk} = testSequence{1, kk}(selFeat, :);
end

testSequence = [];


%% HMM Test
[waveTest.logProb, errors] = mhmm_logprob(waveTest.seq, HMM.priorProb, HMM.trProb,...
                                        HMM.mean, HMM.sigma, HMM.mix);
                                    
%% Evaluate
waveTest.beats = nBeats;
[prob, waveTest.class] = nanmax(waveTest.logProb,[],2);
waveTest.class = waveTest.class(prob > 0);
waveTest.beats = waveTest.beats(prob > 0);
%nb = hist(waveTest.class, 1)


