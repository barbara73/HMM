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
%   - plotAll2Dim.m
%   - plotAll3Dim.m
%   - getType.m
%   - getSum.m
%   - getDeltaT.m
%   - brewermap.m


%% LOAD DATA
clear
close all

addpath(genpath('/Users/barbara/Documents/MATLAB/HMMall'));
pathname = '/Users/barbara/Dropbox/BFH/QT_db/N_patients/';
path = '/Users/barbara/Documents/MATLAB/HMM_QTdb/';

fileID = fopen(strcat(pathname, 'N_PatientsList.txt'));
RECORDS = textscan(fileID, '%s');
fclose(fileID);

indFullyLabelled = {1,2,3,15,16,20,21,5,8};

%for aa = 1 : length(indFullyLabelled)

aa = 1;
kk = indFullyLabelled(aa);
patientName = RECORDS{1, 1}{kk, 1};
load(strcat(path, patientName, '_Data.mat'))
%asynchronous steps - wie am besten speichern??! bestimmte anzahl nehmen?


trainData = {asynchtime, asynchsignal};
labelData = {};


%% Preprocess
% 
% line1 = {'-o','-s','-d','->','-<','-v','-^','-h'};
% line2 = {'--o','--s','--d','-->','--<','--v','--^','--h'};
% marker1 = {'o','s','d','>','<','v','^','h'};
% colorLabel = {'g>','ro','g<'};
% color = {[0 0 1],[1 0 0],[0 1 0],[0 1 1],[1 0 1],[0 0 0.65],[0 0.65 0],[0.65 0 0]};


%% Initialise 
% pathname = '/usr/local/database/qtdb/';


% 
% a = 1:1000;
% b = 1000;

% figure;
% ax(1) = subplot(411)
% plot(tmAsynch(1:b), sigAsynch(1:b), 'Marker', '.')
% title('Asynchronous Sampling')
% ax(2) = subplot(412)
% stem(tmAsynch(1:b), summe(1:b))
% title('Sum of dS over the next 10 steps')
% ax(3) = subplot(413)
% stem(tmAsynch(1:b), du_dt(1:b))
% title('du_dt')
% ax(4) = subplot(414)
% stem(tmAsynch(1:b), count_dt(1:b))
% title('count_dt')
% linkaxes(ax(1:4),'x'); 
% 
% 
% a = 1:1000;
% b = 1000;
% 
% figure;
% ax(1) = subplot(411)
% plot(a, sigAsynch(1:b), 'Marker', '.')
% title('Asynchronous Sampling')
% ax(2) = subplot(412)
% plot(a, summe(1:b))
% title('Sum of dS over the next 10 steps')
% ax(3) = subplot(413)
% plot(a, deltaDS(1:b))
% title('dS(m)*type(m) - dS(m-1)*type(m-1)')
% ax(4) = subplot(414)
% plot(a, dS(1:b))l
% title('dS')
% 
% 
% linkaxes(ax(1:4),'x'); 


% load 19_data_train.mat
% load 19_annotation_train.mat
nbSig = size(asynchtime, 2);
% fs = 1/(t(2)-t(1));

% %% Signal Delineation
% %Init +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% onset = [];
% offset = [];
% peak = [];
% waves = [112,114,116];
% 
% scrsz = get(0,'ScreenSize');
% figure('position',scrsz);
% 
% %Plot signal with existing annotation's
% plot(t,timeSignals(:,1)); hold on;
% if exist('annotation.mat','file');
%     load annotation.mat
%     plot(t(onset(:,1)),timeSignals(onset(:,1),1),colorLabel{1});
%     plot(t(offset(:,1)),timeSignals(offset(:,1),1),colorLabel{3});
%     plot(t(peak(:,1)),timeSignals(peak(:,1),1),colorLabel{2});
% end
% 
% %User Input +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %Button's, key's: 27 = Esc, 30 = Up, 31 = Down, r = 114, p = 112, t = 116
% waitfor(gcf,'CurrentCharacter',30);
% contr = 1;
%     
% while contr ~= 27
%     while contr ~= 31 && contr ~= 27
%         [x,~,button] = ginput(2);
%         if button(1) == 31 || button(1) == 27
%             contr = button(1);
%         elseif button(2) == 31 || button(2) == 27   
%             contr = button(2);
%         elseif button(1) == waves(1) || button(1) == waves(2) || button(1) == waves(3)
%             ind = round(x(2)*fs);
%             if button(2) == 1
%                 plot(t(ind),timeSignals(ind,1),colorLabel{button(2)});
%                 onset = [onset;[ind,find(waves == button(1))]];
%             elseif button(2) == 3
%                 plot(t(ind),timeSignals(ind,1),colorLabel{button(2)});
%                 offset = [offset;[ind,find(waves == button(1))]];
%             elseif button(2) == 2
%                 plot(t(ind),timeSignals(ind,1),colorLabel{button(2)});
%                 peak = [peak;[ind,find(waves == button(1))]];
%             end
%         else
%         end
%     end
%     if contr == 31
%         contr = 0;
%     end
%     while contr == 0
%         waitfor(gcf,'CurrentCharacter',30);
%         contr = 1;
%     end
% end
% 
% %Sort and save ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% [~,ind] = sort(onset(:,1),1);
% onset = onset(ind,:);
% [~,ind] = sort(offset(:,1),1);
% offset = offset(ind,:);
% [~,ind] = sort(peak(:,1),1);
% peak = peak(ind,:);
% 
% save annotation onset offset peak

%% Preprocessing: Wavelet filter
%Wavelet BW & LP Filter +++++++++++++++++++++++++++++++++++++++++++++++++++
levelBW = 9;
levelBeg = 8;
levelEnd = 3;

% %50-Hz Notch-Filter
% wi=2*pi*50/fs;
% p=0.95;
% b=[1 -2*cos(wi) 1];
% a=[1 -2*p*cos(wi) p^2];
waveSig = asynchtime;
for i = 1 : nbSig
    [C,L] = wavedec(asynchtime(:,i),levelBW,'db10');
    BW = wrcoef('a',C,L,'db10',levelBW);      %BW-Estimation = cA(levelBW)              
    D = zeros(1,length(C));                   %Filter High-Frequency noise = cD(1) + cD(2)
    D(sum(L(1:levelBW+2-levelBeg-1))+1:sum(L(1:levelBW+2-levelEnd))) = ...
        C(sum(L(1:levelBW+2-levelBeg-1))+1:sum(L(1:levelBW+2-levelEnd)));
    waveSig(:,i) = waverec(D,L,'db10');
%    waveSig(:,i)= filter(b,a,waveSig(:,i));
end

%1st and 2nd derivative +++++++++++++++++++++++++++++++++++++++++++++++++++
x = waveSig(:,1);
diff1 = x - [x(1);x(1:end-1)];
diff2 = diff1 - [diff1(1);diff1(1:end-1)];

figure;
ax = [];
ax(1) = subplot(2,1,1); hold on;
plot(t,asynchtime(:,1));
plot(t(onset(:,1)),asynchtime(onset(:,1),1),colorLabel{1});
plot(t(offset(:,1)),asynchtime(offset(:,1),1),colorLabel{3});
plot(t(peak(:,1)),asynchtime(peak(:,1),1),colorLabel{2});
ax(2) = subplot(2,1,2); hold on;
plot(t,x);
plot(t,diff1.^2*10,'c');
plot(t,diff2.^2*100,'g');
plot([1,size(x,1)],[0,0],'k');
linkaxes(ax,'x');

%% Feature Selection
%Init +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
lengthVec = 10;
prctVec = 5;
step = lengthVec*prctVec/10;
feat = zeros(size(waveSig,1)/prctVec,nbSig*2);

%Feature generation +++++++++++++++++++++++++++++++++++++++++++++++++++++++
for i = 1:step:size(waveSig,1)-step
        feat(fix(i/prctVec)+1,:) = [sum(waveSig((1:lengthVec)+i-1,:)),...
                  waveSig(i,:)-waveSig(lengthVec+i-1,:)];
end

%Normalization
% feat = feat ./ repmat(repmat(mean(waveSig(peak(peak(:,2) == 2,1),:)-waveSig(peak(peak(:,2) == 1,1),:)),1,2).*...
%                      [repmat(lengthVec,1,4)/2,repmat([1],1,4)],size(feat,1),1);
feat = feat ./ repmat(repmat((max(waveSig,[],1)-min(waveSig,[],1)),1,2).*...
                      [repmat(lengthVec,1,4)/2,repmat([1],1,4)],size(feat,1),1);
% feat = feat ./ repmat(max(feat,[],1),size(feat,1),1);
feat = round(feat*50+50);

%Resample wave indices
onsetFb = [round(onset(:,1)/prctVec),onset(:,2)];
offsetFb = [round(offset(:,1)/prctVec),offset(:,2)];
peakFb = [round(peak(:,1)/prctVec),peak(:,2)];

figure; hold on;
for i = 1:3
    ind = [];
    for j = 1:size(onsetFb,1)
        if onsetFb(j,2) == i && offsetFb(j,2) == i
            ind = [ind;(onsetFb(j,1):offsetFb(j,1))'];
        end
    end
    line = plot(feat(ind,1),feat(ind,5),line1{i});
    set(line,'color',color{i});
end
    ind = [];
    for j = 1:size(offsetFb,1)-1
        if offsetFb(j,2) == 3 && onsetFb(j+1,2) == 1
            ind = [ind;(offsetFb(j,1):onsetFb(j+1,1))'];
        end
    end
    line = plot(feat(ind,1),feat(ind,5),line1{4});
    set(line,'color',color{4});
legend(gca,{'P Wave','R Wave','T Wave','TP Segment'});
     
figure;
ax = [];
ax(1) = subplot(3,1,1); hold on;
sig = decimate(asynchtime(:,1),prctVec);
plot(sig);
plot(round(onset(:,1)/prctVec),sig(round(onset(:,1)/prctVec),1),colorLabel{1});
plot(round(offset(:,1)/prctVec),sig(round(offset(:,1)/prctVec),1),colorLabel{3});
plot(round(peak(:,1)/prctVec),sig(round(peak(:,1)/prctVec),1),colorLabel{2});
ax(2) = subplot(3,1,2); hold on;
plot(feat(:,[1,5]));
plot(onsetFb(:,1),feat(onsetFb(:,1),1),colorLabel{1});
plot(offsetFb(:,1),feat(offsetFb(:,1),1),colorLabel{3});
plot(peakFb(:,1),feat(peakFb(:,1),1),colorLabel{2});
ax(3) = subplot(3,1,3); hold on;
plot(feat(:,[2,6]));
linkaxes(ax,'x');

%% PCA
%Take the most relevant comp ++++++++++++++++++++++++++++++++++++++++++++++
[~,score,latent,~,explained] = pca(feat(:,[1,2,5,6]));
explained
score = round(score*0.9+50);
figure; hold on;
for i = 1:3
    ind = [];
    for j = 1:size(onsetFb,1)
        if onsetFb(j,2) == i && offsetFb(j,2) == i
            ind = [ind;(onsetFb(j,1):offsetFb(j,1))'];
        end
    end
    line = plot(score(ind,1),score(ind,2),line1{i});
    set(line,'color',color{i});
end
    ind = [];
    for j = 1:size(offsetFb,1)-1
        if offsetFb(j,2) == 3 && onsetFb(j+1,2) == 1
            ind = [ind;(offsetFb(j,1):onsetFb(j+1,1))'];
        end
    end
    line = plot(score(ind,1),score(ind,2),line1{4});
    set(line,'color',color{4});
    
figure;
ax = [];
ax(1) = subplot(3,1,1); hold on;
sig = decimate(asynchtime(:,1),prctVec);
plot(sig);
plot(round(onset(:,1)/prctVec),sig(round(onset(:,1)/prctVec),1),colorLabel{1});
plot(round(offset(:,1)/prctVec),sig(round(offset(:,1)/prctVec),1),colorLabel{3});
plot(round(peak(:,1)/prctVec),sig(round(peak(:,1)/prctVec),1),colorLabel{2});
ax(2) = subplot(3,1,2); hold on;
plot(score(:,[1,2]));
plot(onsetFb(:,1),score(onsetFb(:,1),1),colorLabel{1});
plot(offsetFb(:,1),score(offsetFb(:,1),1),colorLabel{3});
plot(peakFb(:,1),score(peakFb(:,1),1),colorLabel{2});
ax(3) = subplot(3,1,3); hold on;
plot(score(:,[3,4]));
linkaxes(ax,'x');

%% HMM Training
%Setup ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Generate sequence in a dataset for 6 different HMMs
HMMs = dataset();
waves = {'PWave','RWave','TWave','PRSegm','RTSegm','TPSegm'};

%Waves (incl one feat samples previous and after)
for i = 1:3
    nbrSeq = min([sum(onsetFb(:,2) == i),sum(offsetFb(:,2) == i)]);
    ind = cell(1,nbrSeq);
    tmp = 1;
    for j = 1:size(onsetFb,1)
        if onsetFb(j,2) == i && offsetFb(j,2) == i
            ind{tmp} = (onsetFb(j,1)-1:offsetFb(j,1)+1)';
            tmp = tmp+1;
        end
    end
    HMMs = [HMMs;dataset(waves(i),nbrSeq,{ind},...
                         'VarNames',{'wave','nbrSeq','seqInd'})];
end

%Segments
for i = 1:3
    if i == 3
        next = 1;
        nbrCorr = 1;
    else
        next = i+1;
        nbrCorr = 0;
    end
    nbrSeq = min([sum(offsetFb(:,2) == i),sum(onsetFb(:,2) == next)])-nbrCorr;
    ind = cell(1,nbrSeq);
    tmp = 1;
    for j = 1:size(offsetFb,1)-nbrCorr
        if offsetFb(j,2) == i && onsetFb(j+1,2) == next
            ind{tmp} = (offsetFb(j,1)+2:onsetFb(j+1,1)-2)';
            tmp = tmp+1;
        end
    end
    HMMs = [HMMs;dataset(waves(i+3),nbrSeq,{ind},...
                         'VarNames',{'wave','nbrSeq','seqInd'})];
end

%Analyze feature vector +++++++++++++++++++++++++++++++++++++++++++++++++++
% figure, hold on;
% for i = 1:size(HMMs,1)
%     seqIndWave = HMMs(i,:).seqInd{1};
%     indAll = [];
%     for j = 1:HMMs.nbrSeq(i)
%         indAll = [indAll;seqIndWave{j}];
%     end
%     ax = subplot(2,3,i);
%     hist(score(indAll,1),50);
%     legend(ax,HMMs(i,:).wave);
% end
% drawnow

%Train HMMs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
waveSel = {'PWave','RWave','TWave'};
Q = [5 4 3 1 1 1];
N = 100;

%Expand dataset
HMMs = [HMMs,dataset({Q','nbrStats'})];
HMMs = [HMMs,dataset({{[],[],[],[],[],[]}','trProb'})];
HMMs = [HMMs,dataset({{[],[],[],[],[],[]}','emitProb'})];
HMMs = [HMMs,dataset({{{},{},{},{},{},{}}','seqStates'})];
HMMs = [HMMs,dataset({{[],[],[],[],[],[]}','seqTests'})];

for i = 1:size(waveSel,2)
    %Extract waves and generate sequences
    sel = strcmp(HMMs.wave,waveSel(i));
    seqIndWave = HMMs(sel,:).seqInd{1};
    seqWave = cell(1,HMMs(sel,:).nbrSeq);
    for j = 1:HMMs(sel,:).nbrSeq
        seqWave{j} = score(seqIndWave{j},[1])';
    end
    
    %Initial parameters
%     trInit = repmat(1/Q,Q,Q);
    trInit = eye(HMMs(sel,:).nbrStats)/2 + circshift(eye(HMMs(sel,:).nbrStats)/2,[-1,0]);
    trInit(end,1) = 0;
    emitInit = repmat(1/N,HMMs(sel,:).nbrStats,N);
    
    %Learn sequences
    [trProb,emitProb] = hmmtrain(seqWave,trInit,emitInit,'Verbose',true);
    HMMs(sel,:).trProb = {trProb};
    HMMs(sel,:).emitProb = {emitProb};
end

%% HMM Visualize
waveSel = {'PWave','RWave','TWave'};

for i = 1:size(waveSel,2)
    %Extract waves and generate sequences
    sel = strcmp(HMMs.wave,waveSel(i));
    seqIndWave = HMMs(sel,:).seqInd{1};
    seqStates = cell(1,HMMs(sel,:).nbrSeq);
    for j = 1:HMMs(sel,:).nbrSeq
        seqStates{j} = hmmviterbi(score(seqIndWave{j},[1]),HMMs(sel,:).trProb{1},HMMs(sel,:).emitProb{1})';
    end  
    HMMs(sel,:).seqStates = {seqStates};
end

figure;
ax = [];
ax(1) = subplot(2,1,1); hold on;
plot(score(:,[1]));
plot(onsetFb(:,1),score(onsetFb(:,1),1),colorLabel{1});
plot(offsetFb(:,1),score(offsetFb(:,1),1),colorLabel{3});
plot(peakFb(:,1),score(peakFb(:,1),1),colorLabel{2});
ax(2) = subplot(2,1,2); hold on;
sig = decimate(asynchtime(:,1),prctVec);
plot(sig(:,[1]),':k');
for i = 1:size(waveSel,2)
    sel = strcmp(HMMs.wave,waveSel(i));
    seqIndWave = HMMs(sel,:).seqInd{1};
    seqStates = HMMs(sel,:).seqStates{1};
    for j = 1:HMMs(sel,:).nbrSeq
        seq = seqIndWave{j};
        states = seqStates{j};
        for ii = 1:size(seq,1)-1
            plot(seq(ii:ii+1),sig(seq(ii:ii+1),[1]),'color',color{states(ii)});
        end
    end
end
linkaxes(ax,'x');

%% HMM Testing
%Init +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load 19_data_train.mat
nbSig = size(asynchtime,2);
fs = 1/(t(2)-t(1));

%Wave Detection +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Setup
levelBW = 9;
levelBeg = 8;
levelEnd = 3;
pointDer = 4;
waveThr = 7e-3;
waveWidth = 0.05;

waveSig = asynchtime;
%Wavelet Bandpass
for i = 1:nbSig
    [C,L] = wavedec(asynchtime(:,i),levelBW,'db10');
    BW = wrcoef('a',C,L,'db10',levelBW);      %BW-Estimation = cA(levelBW)              
    D = zeros(1,length(C));                   %Filter High-Frequency noise = cD(1) + cD(2)
    D(sum(L(1:levelBW+2-levelBeg-1))+1:sum(L(1:levelBW+2-levelEnd))) = ...
        C(sum(L(1:levelBW+2-levelBeg-1))+1:sum(L(1:levelBW+2-levelEnd)));
    waveSig(:,i) = waverec(D,L,'db10');
%    waveSig(:,i)= filter(b,a,waveSig(:,i));
end

%50-Hz Notch-Filter
wi=2*pi*50/fs;
p=0.9;
b=[1 -2*cos(wi) 1];
a=[1 -2*p*cos(wi) p^2];
waveSig= filter(b,a,waveSig);

%Five-point derivative prevents high-frequency noise amplification
diffSig = 1/8*(2*[waveSig(2*pointDer+1:end,:);zeros(2*pointDer,nbSig)] + [waveSig(pointDer+1:end,:);zeros(pointDer,nbSig)] - ...
                [zeros(pointDer,nbSig);waveSig(1:end-pointDer,:)] - 2*[zeros(2*pointDer,nbSig);waveSig(1:end-2*pointDer,:)]);

%Square
sqSig = diffSig.^2;

%Average filter
smSig = reshape(smooth(sqSig,waveWidth*fs),[],nbSig);

%Find global peaks
[peakAmp,peakLoc] = findpeaks(smSig(:,1),'minpeakheight',waveThr,'minpeakdistance',floor(waveWidth*fs));

%Find local maxima
segmLoc = peakLoc(peakLoc-floor(waveWidth/2*fs) >= 1 & peakLoc+floor(waveWidth/2*fs) <= size(waveSig,1)); 
segm = -floor(waveWidth/2*fs)+1:floor(waveWidth/2*fs);
peakSegm = reshape(repmat(segm',1,size(segmLoc,1))+repmat(segmLoc',size(segm,2),1),1,[]);
[~,waveLoc] = max(reshape(abs(waveSig(peakSegm,1)),[],size(segmLoc,1)),[],1);
waveLoc = segmLoc-floor(waveWidth/2*fs) + waveLoc';

figure;
ax = [];
ax(1) = subplot(2,1,1); hold on;
plot(t,asynchtime(:,1));
plot(t,waveSig(:,1),'c');
plot(t(waveLoc),waveSig(waveLoc,1),'ro');
ax(2) = subplot(2,1,2); hold on;
plot(t,sqSig(:,1),'c');
plot(t,smSig(:,1),'g');
plot(t(peakLoc),smSig(peakLoc,1),'ro');
linkaxes(ax,'x');

%Feature Extraction +++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Init
lengthVec = 10;
prctVec = 5;
step = lengthVec*prctVec/10;
feat = zeros(size(waveSig,1)/prctVec,nbSig*2);

%Feature generation
for i = 1:step:size(waveSig,1)-step
    feat(fix(i/prctVec)+1,:) = [sum(waveSig((1:lengthVec)+i-1,:)),...
                  waveSig(i,:)-waveSig(lengthVec+i-1,:)];
end

%Normalization
feat = feat ./ repmat(repmat((max(waveSig,[],1)-min(waveSig,[],1)),1,2).*...
                      [repmat(lengthVec,1,4)/2,repmat([1],1,4)],size(feat,1),1);
feat = round(feat*50+50);

%Decimate peak ind
waveLocFb = round(waveLoc/prctVec);

%PCA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[coeff,score,latent,~,explained] = pca(feat(:,[1,2,5,6]));
explained
featPCA = round(score(:,1)*0.75+50);

%HMM Test +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Setup
waveLength = 25;
waveSel = {'PWave','RWave','TWave'};
logProb = ones(size(waveLocFb,1),size(waveSel,2))*NaN;

%Extract waves around detected peaks
seqIndWave = repmat(-ceil(waveLength/2)+1:floor(waveLength/2),size(waveLocFb,1),1) +...
             repmat(waveLocFb,1,waveLength);

%Test all HHMs for the wave pro p's
for i = 1:size(waveSel,2)
    sel = strcmp(HMMs.wave,waveSel(i));
    for j = 1:size(seqIndWave,1)
        [~,logProb(j,i)] = hmmdecode(featPCA(seqIndWave(j,:))',HMMs(sel,:).trProb{1},HMMs(sel,:).emitProb{1});
    end
end