%% PART 0: Script settings

clear all; close all; clc;

tic

datestring = 'reward';

addpath('.')

subjectInfoFolder = '../../';
eegLocation = '../../results/eeg/';
behaviorFolder = '../../results/behavior/';
saveFolder = '../../results/isc/';
stimLocationAudio = '../../stimuli/';

% read in subject info, see subjects.xlsx file for more info on individual subjects & recording bad channels

subjects = table2struct(readtable([ subjectInfoFolder 'subjects.xlsx']));

[~,z,~] = xlsread([subjectInfoFolder 'subjects.xlsx']);
z = z(2:end,12);

for i = 1:length(subjects)
    if(~isempty(z{i}))
        subjects(i).badChannels = str2num(z{i}); 
    end
end

subjects = subjects(logical([subjects.include]));

EEG = []; % initialize EEG variable to allow concatenation
nComp = 3; % number of components to show
fsRef = 2048;
EEGChannels = 1:64;
EOGChannels = 65:70;
refElectrode = [];
filterCoeff = 3;

% 06/16: changed value of fsDesired from 30 to 64
options.fsDesired = 128;  % desired sampling rate for output (i.e., input into s2e)
options.K = options.fsDesired;  % length of temporal aperture
options.repeat = 1; % number of times to repeat S to handle multiple subjects
options.zscore = 1; % z-score feature
options.highpass = 0; % high-pass feature
% options.skip = (1*options.fsDesired) - 1; % omit this many initial samples
options.skip = 0;
options.numelEEG = 28736;
options.fsRef = fsRef;
options.initialPad = 5 * fsRef;

fsDesired = options.fsDesired;

Nstimuli = 8;

% convert bad channel info to vectors from strings
for i = 1:length(subjects)

  if(isstr(subjects(i).badChannels))
      subjects(i).badChannels = str2num(subjects(i).badChannels);
  end

  if(subjects(i).badChannels == 0)
    subjects(i).badChannels = [];
  end

  load([behaviorFolder subjects(i).behaviorFile]);
  subjects(i).behavior = experimentData;
  subjects(i).conditionOrder = experimentData.conditionOrder;

end

%% PART 1: preprocess stimulus and make stimulus volume

X = []; % initialize stimulus volume

load([stimLocationAudio 'nameOfTheWind003_split_clean.mat']) % load stimulus

for i = Nstimuli:-1:1

    x1 = getEnvelope(audioClean(:,i), fsClean);
    X(:, :, i) = preprocessStimulus(x1, fsClean, options); % preprocesses stimulus, standard preprocessing function in S2R package

end

% semantic dissimilarity processing

% get word embedding
emb = fastTextWordEmbedding;
load('../../stimuli/nameOfTheWindTranscript0003.mat', 'timings');

words = timings;

replacementWords = {'aryen', 'town'; ...
'riverweed', 'plants'; ...
'yeve', 'youve'; ...
'reshi', 'reshi'; ...
'basts', 'bast'; ...
'roah', 'roah'; ...
'waystone', 'lodestone'; ...
'wagoneers', 'drivers'; ...
'dickered', 'haggled'; ...
'hearthfire', 'hearth'; ...
'motherleaf', 'plant'; ...
'tinuë', 'town'; ...
'tarbean', 'town'; ...
'aerueh', 'town'; ...
'chandrian', 'chandrian'; ...
'kvothe', 'kvothe'; ...
'imre', 'town'; ...
'shathered', 'shattered'; ...
'kotes', 'kote'; ...
'nighmane', 'nighmane'; ...
'mhenka', 'mhenka'; ...
'ralien', 'town'; ...
'cealdish', 'cealdish'; ...
'aninn', 'aninn'; ...
'deolan', 'deolan'; ...
'purvis', 'town'; ...
'orrison', 'orrison'; ...
'maneras', 'town'; ...
'groundll', 'ground'};

% replace words not in vocab
for i = 1:height(words)
    if(any(any(strcmp(words{i,'word'}, replacementWords(:,1)))))
        [tmp1] = find(strcmp(words{i,'word'}, replacementWords(:,1)));
        words(i,'word') = replacementWords(tmp1,2);
    end
end

% get embeddings
words.embedding = word2vec(emb,words.word);

% zero out words that arent in vocab and cant be replaced
[ex ix] = find(isnan(words.embedding));
ex = unique(ex);

for i = 1:length(ex)
    words{ex(i),'embedding'} = zeros(1,300);
end

% get which words were presented

load ../../stimuli/splitVidInfo.mat timeSaid

for i = 2:Nstimuli

    timeSaid(:,i:end) = timeSaid(:,i:end) + (options.numelEEG/options.fsDesired);

    timeSaid(find(timeSaid == (options.numelEEG/options.fsDesired))) = 0;

end

for i = 1:length(timeSaid(:))

    [~, tidx] = min(abs(words{:,'start'} - timeSaid(i)));

end

% find stop words

for i = 1:height(words)

    isStopWord(i) = any(strcmp(words{i,3}{1}, stopWords));

end

words.isStopWord = isStopWord';

contentWords = words(~words.isStopWord,:);

% make signal

wordCorrs(1:10)     = 0;
wordCorrsFull(1:10) = 0;

for i = 11:height(contentWords)

    ctemp = corrcoef(mean(contentWords{i-10:i-1,'embedding'},1),contentWords{i,'embedding'});
    if(isnan(ctemp(1,2)))
        % word embeddings
        ctemp(1,2) = 0;
    end
    wordCorrs(i) = 1 - ctemp(1,2);

end

contentWords.semanticDiss = wordCorrs';

for i = 11:height(words)

    ctemp = corrcoef(mean(words{i-10:i-1,'embedding'},1),words{i,'embedding'});
    if(isnan(ctemp(1,2)))
        % word embeddings
        ctemp(1,2) = 0;
    end
    wordCorrsFull(i) = 1 - ctemp(1,2);

end

words.semanticDiss = wordCorrsFull';

semanticSignal     = zeros(length(X)*Nstimuli,1);
semanticSignalFull = zeros(length(X)*Nstimuli,1);
onsetSignal        = zeros(length(X)*Nstimuli,1);
idxSemantic = [1:length(X)*Nstimuli];

for i = 1:height(contentWords)

    [~, ix] = min(abs(idxSemantic - contentWords{i,'start'}*options.fsDesired));
    
    semanticSignal(ix) = contentWords{i,'semanticDiss'};

end

for i = 1:height(words)

    [~, ix] = min(abs(idxSemantic - words{i,'start'}*options.fsDesired));
    
    semanticSignalFull(ix) = words{i,'semanticDiss'};
    onsetSignal(ix) = 1;

end

% split by stimulus
splitIdx = [1:length(X):length(semanticSignal) length(semanticSignal)+1];

for i = 1:Nstimuli

    ss(:,i)  = semanticSignal(splitIdx(i):splitIdx(i+1)-1);
    ssf(:,i) = semanticSignalFull(splitIdx(i):splitIdx(i+1)-1);
    oss(:,i) = onsetSignal(splitIdx(i):splitIdx(i+1)-1);

end

% do toeplitz stuff to the stimulus

semanticSignalProcessed     = [];
semanticSignalProcessedFull = [];
onsetSignalProcessed        = [];

for i = 1:Nstimuli

semanticSignalProcessed(:,:,i)     = cat(2, toeplitz(ss(:,i), [ss(1,i) zeros(1, options.K-2)]), ones(length(X), 1));
semanticSignalProcessedFull(:,:,i) = cat(2, toeplitz(ssf(:,i), [ssf(1,i) zeros(1, options.K-2)]), ones(length(X), 1));
onsetSignalProcessed(:,:,i)        = cat(2, toeplitz(oss(:,i), [oss(1,i) zeros(1, options.K-2)]), ones(length(X), 1));

end

save([saveFolder 'stimVolume_' datestring], 'X', "contentWords", "semanticSignal", 'semanticSignalProcessed', 'semanticSignalProcessedFull',"idxSemantic", 'words', '-v7.3') % save stim volume

% load already processed stim volume
% comment this and uncomment above to process
% load([saveFolder 'stimVolume_reward.mat']);

disp(sprintf('Stimulus processing done \nStarting part 2'))

%% PART 2: Preprocess EEG and make EEG volume

%% preprocess EEG

EEG = [];

% create volume of shape (samples x channels x section x subject)

% load EEG

% for i = 1:length(subjects)
% 
%     load(['../../results/eeg/processed/' subjects(i).eegFile(1:end-4) '.mat'], 'h', 'cutEEG');
% 
%     % initial structure: sample x channel x section
%     cutEEG = double(cutEEG);
%     cutEEG(:,71:end,:) = [];
% 
%     tmp = [];
% 
%     for j = 1:Nstimuli
% 
%         tmp1 = preprocessEEG_v2(cutEEG(:,:,j), EOGChannels, subjects(i).badChannels, ...
%             h.SampleRate, options.fsDesired, options.skip, 5 * h.SampleRate, options.numelEEG);
%         tmp = cat(3, tmp, tmp1);
% 
%     end
% 
%     EEG = cat(4, EEG, tmp);
% 
%     disp(sprintf('done %d', i));
% 
% end
% 
% clear tmp tmp1 cutEEG;
% 
% save('eegVolume', 'EEG', 'subjects', '-v7.3');

% load already processed EEG
% comment this and uncomment above to process instead
load('eegVolume.mat');

% src section
load('eegVolume.mat');

EEG = permute(EEG, [1 3 2 4]);

[T,D,K,N]=size(EEG); % time x trial x channel x subject

EEG = reshape(EEG, [], K, N);

EEG = permute(EEG, [1 3 2]); % swap channels and subjects

EEG = reshape(EEG, [], K);

% reshape stimulus volumes to correct sizes
X = permute(X, [1 3 2]); % rearrange to sample x trial x time lag
X = reshape(X, [], options.K);
X(:, options.K) = 1;

ss = toeplitz(semanticSignal, [semanticSignal(1) zeros(1, options.K-2)]);
ss(:, options.K) = 1;

ssf = toeplitz(semanticSignalFull, [semanticSignalFull(1) zeros(1, options.K-2)]);
ssf(:, options.K) = 1;

oss = toeplitz(onsetSignal, [onsetSignal(1) zeros(1, options.K-2)]);
oss(:, options.K) = 1;

% replicate stimulus for each subject

X = repmat(X, [length(subjects) 1]);
ss = repmat(ss, [length(subjects) 1]);
ssf = repmat(ssf, [length(subjects) 1]);
oss = repmat(oss, [length(subjects) 1]);

%% do SRC

regParam = [10 10];

fprintf(['starting audio src - ' char(datetime('now', 'format', 'hh:mm')) '\n']);
[results.U, results.V, results.H, results.W, results.r, results.p, results.Rxx, results.Ryy]=S2R(X,EEG,regParam);
fprintf(['done audio src - ' char(datetime('now', 'format', 'hh:mm')) '\n']);

fprintf(['starting semantic src - ' char(datetime('now', 'format', 'hh:mm')) '\n']);
[resultsSemantic.U, resultsSemantic.V, resultsSemantic.H, resultsSemantic.W, resultsSemantic.r, resultsSemantic.p, resultsSemantic.Rxx, resultsSemantic.Ryy]=S2R(ss,EEG,regParam);
fprintf(['done semantic src - ' char(datetime('now', 'format', 'hh:mm')) '\n']);

save([saveFolder 'srcResults.mat'], 'results', 'resultsSemantic', 'resultsSemanticFull');

disp('Plotting scalp maps now.')

S2RMap(results.H, results.W, 3, results.r, results.p, options.fsDesired, ...
results.Ryy, results.Rxx, 'BioSemi64.loc');

S2RMap(resultsSemantic.H, resultsSemantic.W, 3, resultsSemantic.r, resultsSemantic.p, options.fsDesired, ...
resultsSemantic.Ryy, resultsSemantic.Rxx, 'BioSemi64.loc');

set(gcf, 'Position', [935 317 699 484]);

saveas(gcf, ['../../images/srcModel_' char(datetime('now', 'format', 'MMddyy-hhmm')) '.png']);

close;

%% make volumes

clear EEG audioClean

load("eegVolume.mat");

load([saveFolder 'stimVolume_reward.mat']);

disp('loaded data');

EEG = permute(EEG, [1 3 2 4]);
X = permute(X, [1 3 2]);

[T,D,K,N]=size(EEG); % time x section x channels x subject

% get sections

eegVolume  = [];
stimVolume = [];

for i = 1:length(subjects)

    c = [];
    d = [];

    for j = 1:4
    
    idx = (subjects(i).conditionOrder == j);
    c(:,:,:,j) = EEG(:, idx, :, i);
    d(:,:,:,j)   = X(:,idx,:);
    
    end

    eegVolume = cat(5, eegVolume, c);
    stimVolume = cat(5, stimVolume, d);

    disp(['done ' num2str(i)]);

end

eegVolume = permute(eegVolume, [1 2 4 3 5]);
stimVolume = permute(stimVolume, [1 2 4 3 5]);

save([saveFolder 'volumeForSRC'], 'eegVolume', 'stimVolume', '-v7.3');

%% do per subject calculation

for i = 1:length(subjects)

    for j = 1:4
    
        tstim = squeeze(stimVolume(:,:,j,:,i));
        tstim = reshape(tstim, [], options.K);

        teeg = squeeze(eegVolume(:,:,j,:,i));
        teeg = reshape(teeg, [], 64); % reshape to time x channels

        U = tstim * results.H;
        V = teeg * results.W;

        nVars=min(size(U,2),size(V,2));
        for n = 1:nVars

            [rhos,pVals]=corrcoef(U(:,n),V(:,n));
            r(n)=rhos(1,2);
            p(n)=pVals(1,2);

        end

        corrs(i,j,:) = r;
        ps(i,j,:) = p;
    
    end

end

srcScores = sum(corrs(:,:,1:3),3);

save([saveFolder 'srcScores'], 'srcScores', '-v7.3');

% anova stuff

noiseLevel  = cat(2, ones(length(srcScores),1), zeros(length(srcScores),1), ones(length(srcScores),1), zeros(length(srcScores),1));
rewardLevel = cat(2, ones(length(srcScores),1), ones(length(srcScores),1), zeros(length(srcScores),1), zeros(length(srcScores),1));
subj        = repmat([1:length(srcScores)]', [1 4]);

[p, t2, stats, terms] = anovan(srcScores(:), {noiseLevel(:), rewardLevel(:), subj(:)}, "random", [3], ...
    varnames={'noise', 'reward', 'subject'}, model='interaction');

%% plot stuff

mywidth = 4;

colors = cat(1, brewermap(NaN,'Set1'), brewermap(NaN,'Set2'), brewermap(NaN,'Set3'), brewermap(NaN,'Dark2'), brewermap(NaN,'Paired'), brewermap(NaN,'Accent'));
colors = cat(2,colors,repmat([0.7],length(colors),1));
colors = repmat(colors, 2, 1);

f = figure;

for i = 1:length(srcScores)

    plot([1 1.5], [srcScores(i,1) srcScores(i,3)], 'color', colors(i,:),'linewidth',mywidth, 'marker', 'o'); hold on;

    plot([2 2.5], [srcScores(i,2) srcScores(i,4)], 'color', colors(i,:),'linewidth',mywidth, 'marker', 'o');

end

means = mean(srcScores);
errors = std(srcScores)/sqrt(length(srcScores));

plot([0.95 1.05],[means(1) means(1)],'-k','LineWidth',3);

plot([1.45 1.55],[means(3) means(3)],'-k','LineWidth',3);

plot([1.95 2.05],[means(2) means(2)],'-k','LineWidth',3);

plot([2.45 2.55],[means(4) means(4)],'-k','LineWidth',3);

plot([1 1.5], means([1 3]), '-k', 'LineWidth', 3);
plot([2 2.5], means([2 4]), '-k', 'LineWidth', 3);

xlim([0.7 2.8]);
xticks([1 1.5 2 2.5])
xticklabels({'-3dB R', '-3dB N', '-9dB R', '-9dB N'});
ylabel('EEG-SRC Score')
title('EEG Results')

[h1 p1] = ttest(srcScores(:,1),srcScores(:,3));
[h2 p2] = ttest(srcScores(:,2), srcScores(:,4));
[h3 p3] = ttest(srcScores(:,1), srcScores(:,2));
[h4 p4] = ttest(srcScores(:,3), srcScores(:,4));

sigstar({[1.25 2.25] [2 2.5]}, [0.001 0.1]);

% ylim([0 1.1])

set(gcf, 'position', [950 497 403 375]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

box off;

saveas(gcf, ['../../images/SRCresults_' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/SRCresults_' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/SRCresults_' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

%% make isc vs. hit pupil fig

saveFolder = '../../results/isc/';

load('behaviorRegressionv2.mat', 'subjects');
load([saveFolder 'srcScores.mat']);
load('missedWordData.mat', 't');
load('extendedPupilResults.mat', 'mpm');

mpm = mpm';
f1 = cat(1,t.f1);

colors = cat(1, brewermap(NaN,'Set1'), brewermap(NaN,'Set2'), brewermap(NaN,'Set3'), brewermap(NaN,'Dark2'), brewermap(NaN,'Paired'), brewermap(NaN,'Accent'));
% colors = cat(2,colors,repmat([0.7],length(colors),1));
% colors = repmat(colors, 2, 1);

hitMeans = [];

for i = 1:length(subjects)

    hitMeans = cat(1, hitMeans, cellfun(@mean, subjects(i).maxHits));

end

f = figure;

% hitMeans for old pupil score, mpm for new extended pupil score (3s pre)
scatter(srcScores(:,1), mpm(:,1), 'Marker', 'square', 'MarkerEdgeColor',colors(1,:),'MarkerEdgeAlpha',0.4); hold on;
scatter(srcScores(:,2), mpm(:,2), 'Marker', 'o', 'MarkerEdgeColor',colors(1,:),'MarkerEdgeAlpha',0.4);
scatter(srcScores(:,3), mpm(:,3), 'Marker', 'square', 'MarkerEdgeColor',colors(2,:),'MarkerEdgeAlpha',0.4);
scatter(srcScores(:,4), mpm(:,4), 'Marker', 'o', 'MarkerEdgeColor',colors(2,:),'MarkerEdgeAlpha',0.4);

hm = mean(mpm);
sm = mean(srcScores);

plot(mean(srcScores(:,1)), hm(1), 'Color', colors(1,:), 'Marker','square', 'MarkerSize', 12, 'MarkerFaceColor',colors(1,:));
plot(mean(srcScores(:,2)), hm(2), 'Color', colors(1,:), 'Marker','o', 'MarkerSize', 12, 'MarkerFaceColor',colors(1,:));
plot(mean(srcScores(:,3)), hm(3), 'Color', colors(2,:), 'Marker','square', 'MarkerSize', 12, 'MarkerFaceColor',colors(2,:));
plot(mean(srcScores(:,4)), hm(4), 'Color', colors(2,:), 'Marker','o', 'MarkerSize', 12, 'MarkerFaceColor',colors(2,:));

hstd = std(mpm)/sqrt(length(subjects));
sstd = std(srcScores)/sqrt(length(subjects));

[e,p] = plotEllipse(sstd(1), hstd(1), sm(1), hm(1), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(2), hstd(2), sm(2), hm(2), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(3), hstd(3), sm(3), hm(3), colors(2,:), 0.3);
[e,p] = plotEllipse(sstd(4), hstd(4), sm(4), hm(4), colors(2,:), 0.3);

xlabel('EEG-SRC');
ylabel('Pre-Target Pupil Size');

set(gcf, 'Position', [997 413 471 414]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

saveas(gcf, ['../../images/SRCvsPupil' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/SRCvsPupil' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/SRCvsPupil' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

% do plot with src vs f1

f = figure;

% hitMeans for old pupil score, mpm for new extended pupil score (3s pre)
scatter(srcScores(:,1), f1(:,1), 'Marker', 'square', 'MarkerEdgeColor',colors(1,:)); hold on;
scatter(srcScores(:,2), f1(:,2), 'Marker', 'o', 'MarkerEdgeColor',colors(1,:));
scatter(srcScores(:,3), f1(:,3), 'Marker', 'square', 'MarkerEdgeColor',colors(2,:));
scatter(srcScores(:,4), f1(:,4), 'Marker', 'o', 'MarkerEdgeColor',colors(2,:));

sm = mean(srcScores);
f1m = mean(f1);

plot(sm(1), f1m(1), 'Color', colors(1,:), 'Marker','square', 'MarkerSize', 12, 'MarkerFaceColor',colors(1,:));
plot(sm(2), f1m(2), 'Color', colors(1,:), 'Marker','o', 'MarkerSize', 12, 'MarkerFaceColor',colors(1,:));
plot(sm(3), f1m(3), 'Color', colors(2,:), 'Marker','square', 'MarkerSize', 12, 'MarkerFaceColor',colors(2,:));
plot(sm(4), f1m(4), 'Color', colors(2,:), 'Marker','o', 'MarkerSize', 12, 'MarkerFaceColor',colors(2,:));

sstd = std(srcScores)/sqrt(length(subjects));
f1std = std(f1)/sqrt(length(subjects));

[e,p] = plotEllipse(sstd(1), f1std(1), sm(1), f1m(1), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(2), f1std(2), sm(2), f1m(2), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(3), f1std(3), sm(3), f1m(3), colors(2,:), 0.3);
[e,p] = plotEllipse(sstd(4), f1std(4), sm(4), f1m(4), colors(2,:), 0.3);

ylabel('F1 Score');
xlabel('EEG-SRC');

% set(gcf, 'Position', [2134 324 666 533]);
set(gcf, 'Position', [997 413 471 414]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

saveas(gcf, ['../../images/F1vsSRC' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/F1vsSRC' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/F1vsSRC' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

% do plot with pupil vs f1

f = figure;

% hitMeans for old pupil score, mpm for new extended pupil score (3s pre)
scatter(mpm(:,1), f1(:,1), 'Marker', 'square', 'MarkerEdgeColor',colors(1,:)); hold on;
scatter(mpm(:,2), f1(:,2), 'Marker', 'o', 'MarkerEdgeColor',colors(1,:));
scatter(mpm(:,3), f1(:,3), 'Marker', 'square', 'MarkerEdgeColor',colors(2,:));
scatter(mpm(:,4), f1(:,4), 'Marker', 'o', 'MarkerEdgeColor',colors(2,:));

hm = mean(mpm);
f1m = mean(f1);

plot(hm(1), f1m(1), 'Color', colors(1,:), 'Marker','square', 'MarkerSize', 12, 'MarkerFaceColor',colors(1,:));
plot(hm(2), f1m(2), 'Color', colors(1,:), 'Marker','o', 'MarkerSize', 12, 'MarkerFaceColor',colors(1,:));
plot(hm(3), f1m(3), 'Color', colors(2,:), 'Marker','square', 'MarkerSize', 12, 'MarkerFaceColor',colors(2,:));
plot(hm(4), f1m(4), 'Color', colors(2,:), 'Marker','o', 'MarkerSize', 12, 'MarkerFaceColor',colors(2,:));

hstd = std(mpm)/sqrt(length(subjects));
f1std = std(f1)/sqrt(length(subjects));

[e,p] = plotEllipse(hstd(1), f1std(1), hm(1), f1m(1), colors(1,:), 0.3);
[e,p] = plotEllipse(hstd(2), f1std(2), hm(2), f1m(2), colors(1,:), 0.3);
[e,p] = plotEllipse(hstd(3), f1std(3), hm(3), f1m(3), colors(2,:), 0.3);
[e,p] = plotEllipse(hstd(4), f1std(4), hm(4), f1m(4), colors(2,:), 0.3);

ylabel('F1 Score');
xlabel('Pre-Target Pupil Size (PTP)');

% set(gcf, 'Position', [2134 324 666 533]);
set(gcf, 'Position', [997 413 471 414]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);


saveas(gcf, ['../../images/F1vsPupil' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/F1vsPupil' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/F1vsPupil' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

%% ALTERNATIVE: DO SUBPLOT WITH THE 3 CHARTS ON IT

f = figure;

markerSize = 7;
q = [];

tiledlayout(2, 7, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile([2 2]);

hm = mean(mpm);
sm = mean(srcScores);

plot(mean(srcScores(:,1)), hm(1), 'Color', colors(1,:), 'Marker','square', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(1,:)); hold on;
plot(mean(srcScores(:,2)), hm(2), 'Color', colors(1,:), 'Marker','o', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(1,:));
plot(mean(srcScores(:,3)), hm(3), 'Color', colors(2,:), 'Marker','square', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(2,:));
plot(mean(srcScores(:,4)), hm(4), 'Color', colors(2,:), 'Marker','o', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(2,:));

hstd = std(mpm)/sqrt(length(subjects));
sstd = std(srcScores)/sqrt(length(subjects));

[e,p] = plotEllipse(sstd(1), hstd(1), sm(1), hm(1), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(2), hstd(2), sm(2), hm(2), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(3), hstd(3), sm(3), hm(3), colors(2,:), 0.3);
[e,p] = plotEllipse(sstd(4), hstd(4), sm(4), hm(4), colors(2,:), 0.3);

xlabel('EEG-SRC');
ylabel('Pre-Target Pupil Size');

set(findall(gcf,'-property','FontSize'),'FontSize',16);

% do plot with src vs f1

nexttile([2 2]);


sm = mean(srcScores);
f1m = mean(f1);

plot(sm(1), f1m(1), 'Color', colors(1,:), 'Marker','square', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(1,:));hold on;
plot(sm(2), f1m(2), 'Color', colors(1,:), 'Marker','o', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(1,:));
plot(sm(3), f1m(3), 'Color', colors(2,:), 'Marker','square', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(2,:));
plot(sm(4), f1m(4), 'Color', colors(2,:), 'Marker','o', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(2,:));

sstd = std(srcScores)/sqrt(length(subjects));
f1std = std(f1)/sqrt(length(subjects));

[e,p] = plotEllipse(sstd(1), f1std(1), sm(1), f1m(1), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(2), f1std(2), sm(2), f1m(2), colors(1,:), 0.3);
[e,p] = plotEllipse(sstd(3), f1std(3), sm(3), f1m(3), colors(2,:), 0.3);
[e,p] = plotEllipse(sstd(4), f1std(4), sm(4), f1m(4), colors(2,:), 0.3);

ylabel('F1 Score');
xlabel('EEG-SRC');

% do plot with pupil vs f1

nexttile([2 2]);

hm = mean(mpm);
f1m = mean(f1);

plot(hm(1), f1m(1), 'Color', colors(1,:), 'Marker','square', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(1,:));hold on;
plot(hm(2), f1m(2), 'Color', colors(1,:), 'Marker','o', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(1,:));
plot(hm(3), f1m(3), 'Color', colors(2,:), 'Marker','square', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(2,:));
plot(hm(4), f1m(4), 'Color', colors(2,:), 'Marker','o', 'MarkerSize', markerSize, 'MarkerFaceColor',colors(2,:));

hstd = std(mpm)/sqrt(length(subjects));
f1std = std(f1)/sqrt(length(subjects));

[e,p] = plotEllipse(hstd(1), f1std(1), hm(1), f1m(1), colors(1,:), 0.3);
[e,p] = plotEllipse(hstd(2), f1std(2), hm(2), f1m(2), colors(1,:), 0.3);
[e,p] = plotEllipse(hstd(3), f1std(3), hm(3), f1m(3), colors(2,:), 0.3);
[e,p] = plotEllipse(hstd(4), f1std(4), hm(4), f1m(4), colors(2,:), 0.3);

ylabel('F1 Score');
xlabel('Pre-Target Pupil Size (PTP)');

set(findall(gcf,'-property','FontSize'),'FontSize',16);

l = legend({'-3dB R', '-9dB R', '-3dB N', '-9dB N'},'Location','layout');
l.Layout.TileSpan = [1 1];
l.Layout.Tile     = 7;

set(gcf, 'Position', [2038 475 1138 362]);

saveas(gcf, ['../../images/F1vsPupilvsSRC_withLine_' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/F1vsPupilvsSRC_withLine_' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/F1vsPupilvsSRC_withLine_' char(datetime('today', 'format', 'MMddyy')) '.fig']);

for i = q; i.Visible = 'off'; end;

saveas(gcf, ['../../images/F1vsPupilvsSRC' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/F1vsPupilvsSRC' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/F1vsPupilvsSRC' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

% do 3d plot with src, f1, and pupil

f1 = cat(1,t.f1);

f = figure;

mpmm = mean(mpm);
sm = mean(srcScores);
fm = mean(f1);

hstd = std(mpm)/sqrt(length(subjects));
sstd = std(srcScores)/sqrt(length(subjects));
fstd = std(f1)/sqrt(length(subjects));

scatter3(srcScores(:,1), mpm(:,1), f1(:,1), 'Marker', 'square', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(1,:)); hold on;
scatter3(srcScores(:,2), mpm(:,2), f1(:,2), 'Marker', 'o'     , 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(1,:));
scatter3(srcScores(:,3), mpm(:,3), f1(:,3), 'Marker', 'square', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(2,:));
scatter3(srcScores(:,4), mpm(:,4), f1(:,4), 'Marker', 'o'     , 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(2,:));

% plot mean and error in 3d ellipse, red for reward, blue for no, darker
% colors are -9 and lighter -3
[x1,y1,z1] = ellipsoid(sm(1), mpmm(1), fm(1), sstd(1), hstd(1), fstd(1));
s1 = surf(x1, y1, z1, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', colors(1,:));

[x2,y2,z2] = ellipsoid(sm(2), mpmm(2), fm(2), sstd(2), hstd(2), fstd(2));
s2 = surf(x2, y2, z2, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#780e0f');

[x3,y3,z3] = ellipsoid(sm(3), mpmm(3), fm(3), sstd(3), hstd(3), fstd(3));
s3 = surf(x3, y3, z3, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', colors(2,:));

[x4,y4,z4] = ellipsoid(sm(4), mpmm(4), fm(4), sstd(4), hstd(4), fstd(4));
s4 = surf(x4, y4, z4, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#1f4666');

xlabel('EEG Score');
ylabel('Pupil Size');
zlabel('F1 Score');

legend({'-3 R', '-9 R', '-3 N', '-9 N'});

% animate fig

nFrames = 30;
view(-61.6800,20.1787);
[az,el] = view;

for i = 1:nFrames

    frame = getframe(gcf);
    img   = frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    if i == 1
        imwrite(img,cmap,'animation.gif','gif','LoopCount',Inf,'DelayTime',0.08);
    else
        imwrite(img,cmap,'animation.gif','gif','WriteMode','append','DelayTime',0.08);
    end

    view(az+i,el);

end

% linear model

subjectVar = repmat([1:20]',[1 4]);
noiseVar   = [zeros(20,1) ones(20,1) zeros(20,1) ones(20,1)];
rewardVar  = [ones(20,2) zeros(20,2)];
conditionVar = repmat([1:4], [20 1]);
varNames = {'subject', 'condition', 'noise', 'reward', 'srcScore', 'pupil', 'f1'};

tbl = table(subjectVar(:), conditionVar(:), nominal(noiseVar(:)), nominal(rewardVar(:)), srcScores(:), mpm(:), f1(:), VariableNames=varNames);

% formula = "f1 ~ 1 + condition + (1 | subject)";
formula = "f1 ~ 1 + noise + reward + srcScore + pupil + (1|subject)";

lme = fitlme(tbl, formula)