%% PART 0: Script settings

clear all; close all; clc;

addpath('.');
addpath('.')
addpath(genpath('../../../analyzeEEG'))
addpath(genpath('../../../StimToResponse'))

tic

datestring = 'pupil';

saveFlag     = 1;
timeSaid     = 1;
separateHits = 1;

mywidth = 4;
colors = cat(1, brewermap(NaN,'Set1'), brewermap(NaN,'Dark2'), brewermap(NaN,'Paired'));

subjectInfoFolder = '../../';
behaviorFolder = '../../results/behavior/';
eegLocation = '../../results/eeg/';
eyeLocation = '../../results/pupil/processed/trimmed/';

saveFolder = '../../results/isc/';

% read in subject info, see subjects.xlsx file for more info on individual subjects & recording bad channels

subjects = table2struct(readtable([ subjectInfoFolder 'subjects.xlsx']));

[~,z,~] = xlsread([subjectInfoFolder 'subjects.xlsx']);
z = z(2:end,12);

for i = 1:length(subjects)
%     if(~(z{i}==0))
    if(~isempty(z{i}))
        subjects(i).badChannels = str2num(z{i}); 
    end
end

subjects = subjects(logical([subjects.include]));
% subjects = subjects(logical([subjects.eyeQuality]));

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

splitSamples = floor(2756525/44100*500);

stimLength   = splitSamples*8;

fsDesired = options.fsDesired;

% load('triggerInfo.mat');
% subjects = t;
% subjects = subjects(logical([subjects.include]));

% convert bad channel info to vectors from strings
for i = 1:length(subjects)

  if(isstr(subjects(i).badChannels))
      subjects(i).badChannels = str2num(subjects(i).badChannels);
  end

  if(subjects(i).badChannels == 0)
    subjects(i).badChannels = [];
  end

  load([behaviorFolder subjects(i).behaviorFile]);
  subjects(i).conditionOrder = experimentData.conditionOrder;

end

%% PART 1: Preprocess pupil and make pupil volume

% load eyes

fsEye = 120;
eyeData = [];

for i = 1:length(subjects)

    load([eyeLocation 'pupil_' num2str(subjects(i).subjectNumber)], 'data', 'breatheData');
    % eyeData = cat(2, eyeData, pupilSignal);
    data2 = squeeze(data(2,:,:));
    data  = squeeze(data(1,:,:));
    subjects(i).pupilData = data;
    subjects(i).rightEyePupil = data2;
    subjects(i).breatheData = fillmissing(breatheData.leftPupil, 'linear');

end

% get baseline
baselines = mean(cat(1,subjects.breatheData),2, 'omitnan');
for i = 1:length(subjects); subjects(i).baseline = baselines(i); end;

% reset length
pupilLength = length(subjects(1).pupilData);

% percent deviation

pdev = [];
ppct = [];

nansignal = [];

filterLen = (fsEye / 5) - 1;

for i = 1:length(subjects)

    t = subjects(i).pupilData;
    q = repmat(t(:,1), [1 length(t)]);
    t = fillmissing(t,'pchip',2);
    t = medfilt1(t ,filterLen, [], 2,'omitnan');
    t = (t-subjects(i).baseline)/subjects(i).baseline*100;
    pdev = cat(3, pdev, t);

    t = subjects(i).pupilData;
    ns = medfilt1(t ,filterLen, [], 2);
    ns(~isnan(ns)) = 0;
    ns(isnan(ns))  = 1;
    nansignal = cat(3,nansignal,ns);
  
    t = subjects(i).pupilData;
    q = repmat(t(:,1), [1 length(t)]);
    t = medfilt1(t ,filterLen, [], 2,'includenan');
    t = fillmissing(t,'linear',2);
    t = (t)/subjects(i).baseline*100;
    ppct = cat(3,ppct, t);

end

save('nansignal.mat', 'nansignal');

% REMOVE THIS TO DO NORMAL PUPIL MEASURES
% ppct = nansignal;

% plot drift figure

x = permute(ppct,[2 1 3]);
y=squeeze(reshape(x,[],1,20));

for i = 1:length(subjects)

    shadedErrorBar([1:length(y)]/fsEye, y(:,i)', repmat(std(y(:,i)),[1 length(y)]), 'lineProps', '-r', 'transparent', 1);

end

% do mean pupil

ppctMean = mean(ppct, 2, 'omitnan');

for i = 1:length(subjects)

    for j = 1:4

        subjects(i).pupilMeans(j) = mean(ppctMean(find(subjects(i).conditionOrder==j),i),'omitnan');

    end

    subjects(i).pupilMeans(find(abs(subjects(i).pupilMeans)>200)) = NaN;

end

f = figure;

for i = 2:length(subjects)

    plot(1, subjects(i).pupilMeans(1), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(1.5, subjects(i).pupilMeans(3), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([1 1.5], [subjects(i).pupilMeans(1) subjects(i).pupilMeans(3)], 'color', colors(i,:),'linewidth',mywidth);

    plot(2, subjects(i).pupilMeans(2), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(2.5, subjects(i).pupilMeans(4), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([2 2.5], [subjects(i).pupilMeans(2) subjects(i).pupilMeans(4)], 'color', colors(i,:),'linewidth',mywidth);

end

xlim([0.7 2.8]);
xticks([1 1.5 2 2.5])
xticklabels({'-3dB R', '-3dB N', '-9dB R', '-9dB N'})
ylabel('Mean Pupil Dilation (% baseline)')
title('Mean Pupil Dilation')

set(gcf, 'position', [950 497 403 375]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

if(saveFlag)
    saveas(gcf, ['../../images/meanPupil_' char(datetime('today', 'format', 'MMddyy')) '.png']);
end

close;

% do hit vs miss analysis

load('missedWordData.mat', 't');

[~,ia,ib] = intersect([subjects.subjectNumber], [t.subjectNumber]);

t = t(ib);

for i = 1:length(t)
    
    subjects(i).behaviorResults   = t(i).behaviorResults;
    subjects(i).adjustedPresses   = t(i).adjustedPresses;
    subjects(i).timeSaid          = t(i).timeSaid       ;
    subjects(i).responseHit       = t(i).responseHit    ;
    subjects(i).missedWords       = t(i).missedWords    ;
    subjects(i).bSignal           = t(i).bSignal        ;
    subjects(i).bSignalPressed    = t(i).bSignalPressed ;
    subjects(i).bSignalPressedHit     = t(i).bSignalPressedHit ;
    subjects(i).bSignalPressedMiss    = t(i).bSignalPressedMiss ;
    
end

% num of s before and after to collect pupil data
pupilWindow = 3; 

for i = 1:length(subjects)

    subjects(i).hitPupil = {};
    subjects(i).missPupil = {};
    
    for j=1:4
        subjects(i).hitPupil{j}  = [];
        subjects(i).missPupil{j} = [];
    end

    for j = 1:length(subjects(i).timeSaid)

        for k = 1:8

            if(subjects(i).timeSaid(j,k) == 0)
                continue;
            end
    
            if(subjects(i).missedWords(j,k) == 0) % ie a hit
    
                ix = floor(subjects(i).timeSaid(j,k)*fsEye);

                if(ix > 26940)
                    continue;
                end
                
                ixs = ix-fsEye*pupilWindow;
                ixe = ix+fsEye*pupilWindow;
                
                if(ixe > 26940)
                    continue;
                end
    
                subjects(i).hitPupil{subjects(i).conditionOrder(k)}(end+1,:) = ppct(k, ixs:ixe, i);
    
            else
    
                ix = floor(subjects(i).timeSaid(j,k)*fsEye);

                if(ix > 26940)
                    continue;
                end

                ixs = ix-fsEye*pupilWindow;
                ixe = ix+fsEye*pupilWindow;
                
                if(ixe > 26940)
                    continue;
                end
    
                subjects(i).missPupil{subjects(i).conditionOrder(k)}(end+1,:) = ppct(k, ixs:ixe, i);

            end

        end

    end

    for j = 1:length(subjects(i).responseHit)

        subjects(i).faPupil{j} = [];

        for k = 1:length(subjects(i).responseHit{j})

            if(subjects(i).responseHit{j}(k)); continue; end;

            ix = floor(subjects(i).adjustedPresses{j}(k)*fsEye);

            if(ix > 26940)
                continue;
            end

            ixs = ix-fsEye*pupilWindow;
            ixe = ix+fsEye*pupilWindow;
            
            if(ixe > 26940 | ixs < 0)
                continue;
            end

            subjects(i).faPupil{j}(end+1,:) = ppct(j,ixs+1:ix,i);

        end

    end

    subjects(i).faPupilMeans = cellfun(@(x) mean(x,2),subjects(i).faPupil,'UniformOutput',false);

end

% get means for period right before hit or miss
for i = 1:length(subjects)
    
    for j = 1:4
       
       subjects(i).beforePresentationMeansHit{j} = mean(subjects(i).hitPupil{j}(:,(pupilWindow*fsEye)-fsEye:pupilWindow*fsEye*1.5),2, 'omitnan')';
       subjects(i).beforePresentationMeansMiss{j} = mean(subjects(i).missPupil{j}(:,(pupilWindow*fsEye)-fsEye:pupilWindow*fsEye*1.5),2, 'omitnan')';
        
    end
    
    subjects(i).rewardVarHit = {ones(1,length(subjects(i).beforePresentationMeansHit{1})), ones(1, length(subjects(i).beforePresentationMeansHit{2})), ...
        zeros(1, length(subjects(i).beforePresentationMeansHit{3})), zeros(1, length(subjects(i).beforePresentationMeansHit{4}))};
    subjects(i).rewardVarMiss = {ones(1,length(subjects(i).beforePresentationMeansMiss{1})), ones(1, length(subjects(i).beforePresentationMeansMiss{2})), ...
        zeros(1, length(subjects(i).beforePresentationMeansMiss{3})), zeros(1, length(subjects(i).beforePresentationMeansMiss{4}))};

    subjects(i).noiseVarHit = {zeros(1,length(subjects(i).beforePresentationMeansHit{1})), ones(1, length(subjects(i).beforePresentationMeansHit{2})), ...
        zeros(1, length(subjects(i).beforePresentationMeansHit{3})), ones(1, length(subjects(i).beforePresentationMeansHit{4}))};
    subjects(i).noiseVarMiss = {zeros(1,length(subjects(i).beforePresentationMeansMiss{1})), ones(1, length(subjects(i).beforePresentationMeansMiss{2})), ...
        zeros(1, length(subjects(i).beforePresentationMeansMiss{3})), ones(1, length(subjects(i).beforePresentationMeansMiss{4}))};
    
end


pupilMeans  = [];
hitOrMiss   = [];
noiseLevel  = [];
rewardLevel = [];
subj = [];

for i = 1:length(subjects)

    pupilMeans  = cat(2, pupilMeans, cat(2,subjects(i).beforePresentationMeansHit{:}), cat(2,subjects(i).beforePresentationMeansMiss{:}));
    hitOrMiss   = cat(2, hitOrMiss, ones(1,length(cat(2,subjects(i).beforePresentationMeansHit{:}))), zeros(1, length(cat(2,subjects(i).beforePresentationMeansMiss{:}))));
    noiseLevel  = cat(2, noiseLevel, cat(2, subjects(i).noiseVarHit{:}), cat(2, subjects(i).noiseVarMiss{:}));
    rewardLevel = cat(2, rewardLevel, cat(2, subjects(i).rewardVarHit{:}), cat(2, subjects(i).rewardVarMiss{:})) ;

    subj        = cat(2, subj, repmat(i, [1 sum([cellfun(@length, subjects(i).beforePresentationMeansHit) cellfun(@length, subjects(i).beforePresentationMeansMiss)])]));

end

% pupilMeans ~ hitOrMiss + noiseLevel + rewardLevel + (1+rewardLevel|subj)

x = table(pupilMeans', hitOrMiss', subj', noiseLevel', rewardLevel', 'VariableNames',{'pupilMeans', 'hitOrMiss', 'subj', 'noiseLevel', 'rewardLevel'});

lme = fitlme(x, 'pupilMeans ~ hitOrMiss + noiseLevel + rewardLevel + (1+rewardLevel|subj)');

[p t2 stats terms] = anovan(pupilMeans, {hitOrMiss subj noiseLevel rewardLevel}, 'random',[2], ...
    'model', 'interaction', 'varnames',{'hitOrMiss', 'subj', 'noiseLevel', 'rewardLevel'}, 'display', 'on');

% do pre button pupil anova

for i = 1:length(subjects)

    for j = 1:length(subjects(i).faPupil)

        if(mod(subjects(i).conditionOrder(j),2)) % if condition # is odd
            nois = 0;
        else
            nois = 1;
        end

        if(subjects(i).conditionOrder(j) < 3) % if condition # is 1 or 2
            rew = 1;
        else
            rew = 0;
        end

        subjects(i).noiseVarFA{j} = repmat(nois, [length(subjects(i).faPupilMeans{j}) 1]);
        subjects(i).rewardVarFA{j} = repmat(rew, [length(subjects(i).faPupilMeans{j}) 1]);

    end

end

% plot hits and misses - aligned to presentation

% plot hits

f = figure; hold on;

for i = 1:length(subjects)

    subplot(2,2,1);

    tem = mean(subjects(i).hitPupil{1},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:721]/fsEye, tem, 'color', colors(i,:)); hold on;

    subplot(2,2,2);

    tem = mean(subjects(i).hitPupil{2},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:721]/fsEye, tem, 'color', colors(i,:)); hold on;

    subplot(2,2,3);

    tem = mean(subjects(i).hitPupil{3},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:721]/fsEye, tem, 'color', colors(i,:)); hold on;

    subplot(2,2,4);

    tem = mean(subjects(i).hitPupil{4},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:721]/fsEye, tem, 'color', colors(i,:)); hold on;

end

for i = 1:4
    
    hold on; subplot(2,2,i);
    xline(361/120);

    xlabel('Time(s)');
    ylabel('Pupil Size (% of baseline)');

    xlim([0 721]/fsEye);
%     ylim([-10 10]);
    ylim([-2 2]);
    xticks([0 120 240 361 481 601 721]/fsEye);
    xticklabels({'-3' '-2' '-1' '0' '1' '2' '3'});
    
    switch i
        case 1
            title('-3dB Reward');
        case 2
            title('-9dB Reward');
        case 3
            title('-3dB No Reward');
        case 4
            title('-9dB No Reward');
    end
    
end

set(gcf, 'position', [950 497 403 375]);
suptitle('Hit Pupil - Spoken Aligned');
saveas(gcf, ['../../images/hitPupil_timeSpokenAligned_' char(datetime('today', 'format', 'MMddyy')) '.png']);

close;

% calculate deflection after word spoken

for i = 1:length(subjects)

    for j = 1:4
        subjects(i).maxHits{j} = max(subjects(i).hitPupil{j}(:,361:501),[],2);
    end

end

maxHits = cell(1,4);
subjVar = cell(1,4);
noiseVar = cell(1,4);
rewardVar = cell(1,4);

for i = 1:length(subjects)

    for j = 1:4

        % before word
        subjects(i).hitPupil{j} = subjects(i).hitPupil{j}(:,1:3*fsEye);
        subjects(i).missPupil{j} = subjects(i).missPupil{j}(:,1:3*fsEye);

        % after word
%         subjects(i).hitPupil{j} = subjects(i).hitPupil{j}(:,3*fsEye:end);
%         subjects(i).missPupil{j} = subjects(i).missPupil{j}(:,3*fsEye:end);

    end

end

for i = 1:length(subjects)

    subjects(i).meanPupilHit = cellfun(@(x) mean(x,2), subjects(i).hitPupil, 'UniformOutput',false);
    subjects(i).meanPupilMiss = cellfun(@(x) mean(x,2), subjects(i).missPupil, 'UniformOutput',false);

    for j = 1:4

        maxHits{j} = cat(1,maxHits{j}, subjects(i).maxHits{j});
        subjVar{j} = cat(1,subjVar{j}, repmat([i], [length(subjects(i).maxHits{j}) 1]));

    end
end

noiseVar{1} = zeros(length(maxHits{1}),1);
noiseVar{2} = ones(length(maxHits{2}),1);
noiseVar{3} = zeros(length(maxHits{3}),1);
noiseVar{4} = ones(length(maxHits{4}),1);

rewardVar{1} = ones(length(maxHits{1}),1);
rewardVar{2} = ones(length(maxHits{2}),1);
rewardVar{3} = zeros(length(maxHits{3}),1);
rewardVar{4} = zeros(length(maxHits{4}),1);

subjVar = cat(1,subjVar{:});
noiseVar = cat(1,noiseVar{:});
rewardVar = cat(1, rewardVar{:});
maxHits = cat(1, maxHits{:});

[p t stats] = anovan(maxHits, {noiseVar, rewardVar, subjVar}, 'random', [3], ...
    'model', 'interaction', 'varnames', {'noise', 'reward', 'subj'}, 'display', 'on');

% do mean across extended pupil - +/- 3s from word presentation

mp = zeros(4, length(subjects), 2); % condition x subject x hit/miss
s = zeros(4, length(subjects), 2);  % 1 for hit, 2 for miss
n = zeros(4, length(subjects), 2);
r = zeros(4, length(subjects), 2);

for i = 1:length(subjects)

    subjects(i).meanPupilHit  = cellfun(@(x) mean(x,2), subjects(i).hitPupil, 'UniformOutput',false);
    subjects(i).meanPupilMiss = cellfun(@(x) mean(x,2), subjects(i).missPupil, 'UniformOutput',false);

    mp(:,i,1) = cellfun(@mean, subjects(i).meanPupilHit);
    mp(:,i,2) = cellfun(@mean, subjects(i).meanPupilMiss);
    s(:,i,:)  = i;

end

n([2 4],:,:) = 1;
r([1 2],:,:) = 1;

[p t stats] = anovan(mp(:), {n(:), r(:), s(:)}, 'random', [3], ...
    'model', 'interaction', 'varnames', {'noise', 'reward', 'subj'});

% plot extended pupil mean

colors = cat(1, brewermap(NaN,'Set1'), brewermap(NaN,'Set2'), brewermap(NaN,'Set3'), brewermap(NaN,'Dark2'), brewermap(NaN,'Paired'), brewermap(NaN,'Accent'));
colors = cat(2,colors,repmat([0.7],length(colors),1));
colors = repmat(colors, 2, 1);

f = figure; hold on;

mpm = mean(mp,3); % take average of hits and misses for each sub/cond

for i = 1:length(subjects)

    plot([1 1.5], [mpm(1,i) mpm(3,i)], 'color', colors(i,:),'linewidth',mywidth, 'marker', 'o');
    plot([2 2.5], [mpm(2,i) mpm(4,i)], 'color', colors(i,:),'linewidth',mywidth, 'marker', 'o');

end

means = mean(mpm');
errors = std(mpm')/sqrt(length(mpm'));

plot([0.95 1.05],[means(1) means(1)],'-k','LineWidth',3);
% patch([0.95 1.05 1.05 0.95], [means(1)-errors(1) means(1)-errors(1), means(1)+errors(1) means(1)+errors(1)], [0 0 0], 'facealpha',0.3)

plot([1.45 1.55],[means(3) means(3)],'-k','LineWidth',3);
% patch([1.45 1.55 1.55 1.45], [means(3)-errors(3) means(3)-errors(3), means(3)+errors(3) means(3)+errors(3)], [0 0 0], 'facealpha',0.3)

plot([1.95 2.05],[means(2) means(2)],'-k','LineWidth',3);
% patch([1.95 2.05 2.05 1.95], [means(2)-errors(2) means(2)-errors(2), means(2)+errors(2) means(2)+errors(2)], [0 0 0], 'facealpha',0.3)

plot([2.45 2.55],[means(4) means(4)],'-k','LineWidth',3);
% patch([2.45 2.55 2.55 2.45], [means(4)-errors(4) means(4)-errors(4), means(4)+errors(4) means(4)+errors(4)], [0 0 0], 'facealpha',0.3)

plot([1 1.5], means([1 3]), '-k', 'LineWidth', 3);
plot([2 2.5], means([2 4]), '-k', 'LineWidth', 3);

xlim([0.7 2.8]);
xticks([1 1.5 2 2.5])
xticklabels({'-3dB R', '-3dB N', '-9dB R', '-9dB N'});
ylabel('Pupil Size (% of baseline)')
title('Pupil Results')

% [h1 p1] = ttest(mpm(1,:), mpm(3,:));
% [h2 p2] = ttest(mpm(2,:), mpm(4,:));
% [h3 p3] = ttest(mpm(1,:), mpm(2,:));
% [h4 p4] = ttest(mpm(3,:), mpm(4,:));

[COMPARISON,MEANS,H,GNAMES] = multcompare(stats,'dimension',[1 2],Display='off');

p1 = COMPARISON(2,6);
p2 = COMPARISON(5,6);
p3 = COMPARISON(6,6);
p4 = COMPARISON(1,6);

% sigstar({[1 1.5] [2 2.5] [1 2] [1.5 2.5]}, [p1 p2 p3 p4]);
sigstar({[1.25 2.25] [2 2.5]}, [0.5 p2])

set(gcf, 'position', [950 497 403 375]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

box off;

saveas(gcf, ['../../images/pupilBySubjectCondition_' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/pupilBySubjectCondition_' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/pupilBySubjectCondition_' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

% by trial calculation
% meanPupil = cell(1,4);
% subjVar   = cell(1,4);
% noiseVar  = cell(1,4);
% rewardVar = cell(1,4);
% 
% for i = 1:length(subjects)
% 
%     subjects(i).meanPupilHit  = cellfun(@(x) mean(x,2), subjects(i).hitPupil, 'UniformOutput',false);
%     subjects(i).meanPupilMiss = cellfun(@(x) mean(x,2), subjects(i).missPupil, 'UniformOutput',false);
% 
%     for j = 1:4
% 
%         meanPupil{j} = cat(1, meanPupil{j}, subjects(i).meanPupilHit{j}, subjects(i).meanPupilMiss{j});
%         subjVar{j} = cat(1,subjVar{j}, repmat([i], [length(subjects(i).meanPupilHit{j})+length(subjects(i).meanPupilMiss{j}) 1]));
% 
%     end
% end
% 
% noiseVar{1} = zeros(length(meanPupil{1}),1);
% noiseVar{2} = ones(length(meanPupil{2}),1);
% noiseVar{3} = zeros(length(meanPupil{3}),1);
% noiseVar{4} = ones(length(meanPupil{4}),1);
% 
% rewardVar{1} = ones(length(meanPupil{1}),1);
% rewardVar{2} = ones(length(meanPupil{2}),1);
% rewardVar{3} = zeros(length(meanPupil{3}),1);
% rewardVar{4} = zeros(length(meanPupil{4}),1);
% 
% subjVar = cat(1,subjVar{:});
% noiseVar = cat(1,noiseVar{:});
% rewardVar = cat(1, rewardVar{:});
% meanPupil = cat(1, meanPupil{:});
% 
% [p t stats] = anovan(meanPupil, {noiseVar, rewardVar, subjVar}, 'random', [3], ...
%     'model', 'interaction', 'varnames', {'noise', 'reward', 'subj'});

% plot max hit per subject and condition

colors = cat(1, brewermap(NaN,'Set1'), brewermap(NaN,'Set2'), brewermap(NaN,'Set3'), brewermap(NaN,'Dark2'), brewermap(NaN,'Paired'), brewermap(NaN,'Accent'));
colors = cat(2,colors,repmat([0.7],length(colors),1));
colors = repmat(colors, 2, 1);

f = figure; hold on;

for i = 1:length(subjects)

    mh = cellfun(@mean, subjects(i).maxHits);

    plot([1 1.5], [mh(1) mh(3)], 'color', colors(i,:),'linewidth',mywidth, 'marker', 'o');
    plot([2 2.5], [mh(2) mh(4)], 'color', colors(i,:),'linewidth',mywidth, 'marker', 'o');

end

xlim([0.7 2.8]);
xticks([1 1.5 2 2.5])
xticklabels({'-3dB R', '-3dB N', '-9dB R', '-9dB N'});
ylabel('Pupil Size (% of baseline)')
title('Pupil Results')

mh = cat(1,subjects.maxHits);
mh = cellfun(@mean, mh);

[h1 p1] = ttest(mh(:,1), mh(:,3));
[h2 p2] = ttest(mh(:,2), mh(:,4));
[h3 p3] = ttest(mh(:,1), mh(:,2));
[h4 p4] = ttest(mh(:,3), mh(:,4));

sigstar({[1 1.5] [2 2.5] [1 2] [1.5 2.5]}, [p1 p2 p3 p4]);

set(gcf, 'position', [950 497 403 375]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

box off;

saveas(gcf, ['../../images/hitPupilAcrossSubjects' char(datetime('today', 'format', 'MMddyy')) '.png']);
saveas(gcf, ['../../images/hitPupilAcrossSubjects' char(datetime('today', 'format', 'MMddyy')) '.epsc']);
saveas(gcf, ['../../images/hitPupilAcrossSubjects' char(datetime('today', 'format', 'MMddyy')) '.fig']);

close;

save('behaviorRegressionv2.mat', 'subjects', '-v7.3');

% plot misses

f = figure; hold on;

for i = 1:length(subjects)

    subplot(2,2,1);

    tem = mean(subjects(i).missPupil{1},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:length(tem)]/fsEye, tem, 'color', colors(i,:)); hold on;

    subplot(2,2,2);

    tem = mean(subjects(i).missPupil{2},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:length(tem)]/fsEye, tem, 'color', colors(i,:)); hold on;

    subplot(2,2,3);

    tem = mean(subjects(i).missPupil{3},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:length(tem)]/fsEye, tem, 'color', colors(i,:)); hold on;

    subplot(2,2,4);

    tem = mean(subjects(i).missPupil{4},1, 'omitnan');
    tem = tem - mean(tem);
    plot([1:length(tem)]/fsEye, tem, 'color', colors(i,:)); hold on;

end

for i = 1:4
    
    hold on; subplot(2,2,i);
    xline(361/120);

    xlabel('Time(s)');
    ylabel('Pupil Size (% of baseline)');

    xlim([0 721]/fsEye);
%     ylim([-10 10]);
    ylim([-2 2]);
    xticks([0 120 240 361 481 601 721]/fsEye);
    xticklabels({'-3' '-2' '-1' '0' '1' '2' '3'});
    
    switch i
        case 1
            title('-3dB Reward');
        case 2
            title('-9dB Reward');
        case 3
            title('-3dB No Reward');
        case 4
            title('-9dB No Reward');
    end
    
end

set(gcf, 'position', [2366 87 1243 798]);
suptitle('Miss Pupil - Spoken Aligned');
saveas(gcf, ['../../images/missPupil_timeSpokenAligned_' char(datetime('today', 'format', 'MMddyy')) '.png']);

close;

% extract pupil aligned to button press

for i = 1:length(subjects)

    sig = [1:length(subjects(i).pupilData)]/fsEye;
    subjects(i).responseAlignedPupil = cell(1,8);

    for j = 1:8

        subjects(i).responseAlignedPupil{j} = zeros(length(subjects(i).adjustedPresses{j}), (6*fsEye)+1);

        for k = 1:length(subjects(i).adjustedPresses{j})

            [x ix] = min(abs(sig - subjects(i).adjustedPresses{j}(k)));

            ixs = ix - fsEye*3;
            ixe = ix + fsEye*3;

            if(ixs < 1 | ixe > length(ppct))
                subjects(i).responseAlignedPupil{j}(k,:) = NaN;
                disp(fprintf('subject %d, section %d, response %d could not be extracted', i, j, k));
                continue;
            end

            subjects(i).responseAlignedPupil{j}(k,:) = ppct(j, ixs:ixe, i);

        end

    end

end

% plot pupil aligned to button press - hits only

f = figure; hold on;

for i = 1:length(subjects)

    tem = cat(1, subjects(i).responseAlignedPupil{:});
    idx = logical(cat(2, subjects(i).responseHit{:}));

    tem = tem(idx, :);

    tem = mean(tem, 1, 'omitnan');
    tem = tem - mean(tem);

    plot([1:721]/fsEye, tem, 'color', colors(i,:));

end

xline(361/120);

xlabel('Time(s)');
ylabel('Pupil Size (% of baseline)');

xlim([0 721]/fsEye);
% ylim([-10 10]);
ylim([-2 2]);
xticks([0 120 240 361 481 601 721]/fsEye);
xticklabels({'-3' '-2' '-1' '0' '1' '2' '3'});

title('Pupil Size Aligned to Button Press - Hits');

set(gcf, 'position', [2214 281 758 594]);
saveas(gcf, ['../../images/hitPupil_buttonAligned_' char(datetime('today', 'format', 'MMddyy')) '.png']);

close;

% plot pupil aligned to button press - false alarm only

f = figure; hold on;

for i = 1:length(subjects)

    tem = cat(1, subjects(i).responseAlignedPupil{:});
    idx = ~logical(cat(2, subjects(i).responseHit{:}));

    tem = tem(idx, :);

    tem = mean(tem, 1, 'omitnan');
    tem = tem - mean(tem);

    plot([1:721]/fsEye, tem, 'color', colors(i,:));

    disp(fprintf('subject %d has %d false alarms', subjects(i).subjectNumber, length(find(idx))));

end

xline(361/120);

xlabel('Time(s)');
ylabel('Pupil Size (% of baseline)');

xlim([0 721]/fsEye);
% ylim([-10 10]);
ylim([-2 2]);
xticks([0 120 240 361 481 601 721]/fsEye);
xticklabels({'-3' '-2' '-1' '0' '1' '2' '3'});

title('Pupil Size Aligned to Button Press - False Alarm');

set(gcf, 'position', [2214 281 758 594]);
saveas(gcf, ['../../images/falseAlarmPupil_buttonAligned_' char(datetime('today', 'format', 'MMddyy')) '.png']);

close;

% hits minus false alarm

f = figure; hold on;

for i = 1:length(subjects)

    tem = cat(1, subjects(i).responseAlignedPupil{:});
    idx = ~logical(cat(2, subjects(i).responseHit{:}));

    tem = tem(idx, :);

    tem = mean(tem, 1, 'omitnan');
    tem = tem - mean(tem);

    tem2 = cat(1, subjects(i).hitPupil{:});
    tem2 = mean(tem2, 1, 'omitnan');
    tem2 = tem2 - mean(tem2);

    plot([1:length(tem)]/fsEye, tem2 - tem, 'color', colors(i,:));

end

xline(361/120);

xlabel('Time(s)');
ylabel('Pupil Size (% of baseline)');

xlim([0 721]/fsEye);
ylim([-10 10]);
xticks([0 120 240 361 481 601 721]/fsEye);
xticklabels({'-3' '-2' '-1' '0' '1' '2' '3'});

title('Pupil Size Aligned to Button Press - Hit minus FA');

set(gcf, 'position', [2214 281 758 594]);
saveas(gcf, ['../../images/subtractionPupil_hitMinusFalseAlarm_' char(datetime('today', 'format', 'MMddyy')) '.png']);

close;

% get response aligned by condition

for i = 1:length(subjects)

    subjects(i).responseAlignedPupilByCondition = [];

    for j = 1:4

        idx = find(subjects(i).conditionOrder == j);

        subjects(i).responseAlignedPupilByCondition{j} = cat(1, subjects(i).responseAlignedPupil{idx});

    end

end

for i = 1:length(subjects)

    t1 = []; t2 = [];

    for j = 1:4

        t2(j,:) = mean(subjects(i).responseAlignedPupilByCondition{j}, 1, 'omitnan');

        t1(j,:) = mean(subjects(i).hitPupil{j}, 1, 'omitnan');

    end

    subjects(i).pupilSubtraction = t1 - t2;

    t = mean(subjects(i).pupilSubtraction, 2, 'omitnan');
    subjects(i).pupilSubtraction = subjects(i).pupilSubtraction - repmat(t, [1 721]);

end

% plot pupil subtraction

f = figure;

for i = 1:length(subjects)

    for j = 1:4

        subplot(2,2,j); hold on;

        plot([1:721]/fsEye, subjects(i).pupilSubtraction(j,:), 'color', colors(i,:)); hold on;


    end

end

for i = 1:4
    
    hold on; subplot(2,2,i);
    xline(361/120);

    xlabel('Time(s)');
    ylabel('Pupil Size (% of baseline)');

    xlim([0 721]/fsEye);
    xticks([0 120 240 361 481 601 721]/fsEye);
    xticklabels({'-3' '-2' '-1' '0' '1' '2' '3'});
    
    switch i
        case 1
            title('-3dB Reward');
        case 2
            title('-9dB Reward');
        case 3
            title('-3dB No Reward');
        case 4
            title('-9dB No Reward');
    end
    
end

set(gcf, 'position', [2366 87 1243 798]);
suptitle('Subtraction of Hit Pupil & Button Aligned Pupil');
saveas(gcf, ['../../images/subtractionOfHitPupilAndButtonAligned_' char(datetime('today', 'format', 'MMddyy')) '.png']);

close;

% calculate src

% make toeplitz

% K = 5*fsEye;

fsRegress = 60;
K = 10 * fsRegress;

for i = 1:length(subjects)

    c1 = [subjects(i).bSignal(1,1) zeros(1,K-1)];
    tmp1 = toeplitz(subjects(i).bSignal(:,1), c1);

    c2 = [subjects(i).bSignalPressedHit(1,1) zeros(1,K-1)];
    tmp2 = toeplitz(subjects(i).bSignalPressedHit(:,1), c1);

    c3 = [subjects(i).bSignalPressedMiss(1,1) zeros(1,K-1)];
    tmp3 = toeplitz(subjects(i).bSignalPressedMiss(:,1), c1);

    for j = 2:8

        c1 = [subjects(i).bSignal(1,j) zeros(1,K-1)];
        c2 = [subjects(i).bSignalPressedHit(1,j) zeros(1,K-1)];
        c3 = [subjects(i).bSignalPressedMiss(1,j) zeros(1,K-1)];

        tmp1 = cat(3, tmp1, toeplitz(subjects(i).bSignal(:,j), c1));
        tmp2 = cat(3, tmp2, toeplitz(subjects(i).bSignalPressedHit(:,j), c2));
        tmp3 = cat(3, tmp3, toeplitz(subjects(i).bSignalPressedMiss(:,j), c2));

    end

    subjects(i).bSignaltplz = tmp1;
    subjects(i).bSignalPressedHittplz = tmp2;
    subjects(i).bSignalPressedMisstplz = tmp3;

end

clear tmp1 tmp2 tmp3;

% concatenate behavior signals

bSignal = [];
bSignalPressedHit = [];
bSignalPressedMiss = [];

for i = 1:length(subjects)

    subjects(i).bSignaltplz = permute(subjects(i).bSignaltplz, [1 3 2]);
    subjects(i).bSignalPressedHittplz = permute(subjects(i).bSignalPressedHittplz, [1 3 2]);
    subjects(i).bSignalPressedMisstplz = permute(subjects(i).bSignalPressedMisstplz, [1 3 2]);

    bSignal = cat(1, bSignal, reshape(subjects(i).bSignaltplz, [], K));
    bSignalPressedHit = cat(1, bSignalPressedHit, reshape(subjects(i).bSignalPressedHittplz, [], K));
    bSignalPressedMiss = cat(1, bSignalPressedMiss, reshape(subjects(i).bSignalPressedMisstplz, [], K));

end

% src

bSignal(find(bSignal == -1)) = 1;

bSignal = cat(2, bSignal, bSignalPressedHit, bSignalPressedMiss);

bSignal = cat(2, bSignal, ones(length(bSignal),1));

clear bSignalPressedHit bSignalPressedMiss

y = permute(ppct, [2 1 3]);
y = downsample(y,2);

H1 = bSignal\y(:);

[S T N] = size(y);
zeropad = zeros(fsRegress*5, T, N); 
y = cat(1, zeropad, y);
y = y(1:S, :, :);

H2 = bSignal\y(:);

y = permute(ppct, [2 1 3]);
y = downsample(y,2);
zeropad = zeros(fsRegress*5, T, N); 
y = cat(1, zeropad, y);
y = y(1:S, :, :);

for i = 1:length(subjects)

    for j = 1:4

        [~, idx] = find(subjects(i).conditionOrder == j);

        % bSigPerSubject = cat(2, reshape(subjects(i).bSignaltplz, [], K), ...
        %     reshape(subjects(i).bSignalPressedHittplz, [], K), reshape(subjects(i).bSignalPressedMisstplz, [], K));

%         bSigPerSubject = cat(2, reshape(subjects(i).bSignaltplz(:, idx, :), [], K), ...
%             reshape(subjects(i).bSignalPressedHittplz(:, idx, :), [], K));

        bSigPerSubject = cat(2, reshape(subjects(i).bSignaltplz(:, idx, :), [], K), ...
            reshape(subjects(i).bSignalPressedHittplz(:, idx, :), [], K) + reshape(subjects(i).bSignalPressedMisstplz(:, idx, :), [], K), ...
            reshape(subjects(i).bSignalPressedMisstplz(:, idx, :), [], K));

        bSigPerSubject = cat(2, bSigPerSubject, ones(length(bSigPerSubject),1));

        yPerSub = y(:, idx, i);

        subjects(i).H(:,j) = bSigPerSubject\yPerSub(:);

        disp(['done sub ' num2str(subjects(i).subjectNumber) ' condition ' num2str(j)]);

    end

end

save('behaviorRegressionPerSubject', 'H1', 'H2', 'subjects', '-v7.3');

keyboard;

% [U, V, H, W, r, p, Rxx, Ryy] = S2R(bSignal,y(:), [1 1]);

% save('srcBehavioralData', 'U', 'V', 'H', 'W', 'r', 'p', 'Rxx', 'Ryy', '-v7.3');

% EEG section

% do isc

ppct = permute(ppct, [2 1 3]);
x = reshape(ppct, [], length(subjects));
x(isnan(x)) = 100;
X = [];
X(:,1,:) = x;

clear x;

[isc, sub, sec, w, a] = isceeg(X, fsEye);

% calculate isc

x = [];
x(:,1,:) = ppct;
[isc, sub, sec, w, a] = isceeg(x, fsEye);

h = cat(2,squeeze(hardSections(:,1,:)),squeeze(hardSections(:,2,:)),squeeze(hardSections(:,3,:)),squeeze(hardSections(:,4,:)));
e = cat(2,squeeze(easySections(:,1,:)),squeeze(easySections(:,2,:)),squeeze(easySections(:,3,:)),squeeze(easySections(:,4,:)));

hs(:,1,:) = h;
es(:,1,:) = e;

hs = permute(hs,[3 2 1]);
es = permute(es,[3 2 1]);

[T,D,N] = size(hs);

Rij_hard = permute(reshape(cov(hs(:,:)),[D N  D N]),[1 3 2 4]);

Rij_easy = permute(reshape(cov(es(:,:)),[D N  D N]),[1 3 2 4]);

% compute hard
for i=1:N
    Rw=0; for j=1:N, if i~=j, Rw = Rw+1/(N-1)*(Rij_hard(:,:,i,i)+Rij_hard(:,:,j,j)); end; end
    Rb=0; for j=1:N, if i~=j, Rb = Rb+1/(N-1)*(Rij_hard(:,:,i,j)+Rij_hard(:,:,j,i)); end; end
    ISCpersubject_hard(:,i) = diag(w'*Rb*w)./diag(w'*Rw*w);
end

% compute easy
for i=1:N
    Rw=0; for j=1:N, if i~=j, Rw = Rw+1/(N-1)*(Rij_easy(:,:,i,i)+Rij_easy(:,:,j,j)); end; end
    Rb=0; for j=1:N, if i~=j, Rb = Rb+1/(N-1)*(Rij_easy(:,:,i,j)+Rij_easy(:,:,j,i)); end; end
    ISCpersubject_easy(:,i) = diag(w'*Rb*w)./diag(w'*Rw*w);
end

isce = ISCpersubject_easy;
isch = ISCpersubject_hard;

f = figure;

for i = 1:length(isch)

    plot(1.5, isch(i), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(1, isce(i), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([1 1.5], [isce(i) isch(i)], 'color', colors(i,:),'linewidth',mywidth);

end

xlim([0.7 1.8]);
xticks([1 1.5])
xticklabels({'-3dB', '-9dB'})
ylabel('ISC Score')
title('ISC Results (Pupil)')

[h p] = ttest(isch,isce);

sigstar({[1 1.5]}, p);

% ylim([0 1.1])

set(gcf, 'position', [680 484 492 494]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

if(saveFlag)
    saveas(gcf, ['../../images/pupilISC_' char(datetime('today', 'format', 'MMddyy')) '.png']);
end

% display means per subject/section

f = figure;

for i = 1:length(easyMean)
    
    for j = 1:4

        plot(1.5, hardMean(i,j), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
        plot(1, easyMean(i,j), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
        plot([1 1.5], [easyMean(i,j) hardMean(i,j)], 'color', colors(i,:),'linewidth',mywidth);
    
    end

end

xlim([0.7 1.8]);
xticks([1 1.5])
xticklabels({'-3dB', '-9dB'})
ylabel('Pupil Dilation (% baseline)')
title('Mean pupil per subject/section')

[h p] = ttest(hardMean(:),easyMean(:));

sigstar({[1 1.5]}, p);

% ylim([0 1.1])

set(gcf, 'position', [680 484 492 494]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);

if(saveFlag)
    saveas(gcf, ['../../images/meanPupilBySubjectSection_' char(datetime('today', 'format', 'MMddyy')) '.png']);
end

close;

% display 1,2,3 min words for subject 13
displayTimes = [60 120 180];

% strace = pdev(:,3);
strace = mean(pdev,2);

close all;

for i = 1:length(displayTimes)
    
    figure;
    
    b = (displayTimes(i)-5)*fsEye;
    e = (displayTimes(i)+5)*fsEye;
    e = e-1;
    
    plot([1:length(strace(b:e))]/500,strace(b:e)); hold on;
    scatter(5, strace(b+2500),60,'k','filled');
    
    xticks([0 5 10]);
    xticklabels({'-5','0','+5'});
    
    xlabel('Seconds (s)');
    ylabel('Pupil Dilation');
    
    if(i==1)
        title('Word: road');
        if(saveFlag); saveas(gcf, ['../../images/wordTrace1_' char(datetime('today', 'format', 'MMddyy')) '.png']); end;
    elseif(i==2)
        title('Word: inside');
        if(saveFlag); saveas(gcf, ['../../images/wordTrace2_' char(datetime('today', 'format', 'MMddyy')) '.png']); end;
    else
        title('Word: spoke');
        if(saveFlag); saveas(gcf, ['../../images/wordTrace3_' char(datetime('today', 'format', 'MMddyy')) '.png']); end;
    end
    
end

close all

% calculate mean values across words using window

subjectData = t;

windowSize = 5;

p = ppct;

wordsplits = [1:fsEye*5:length(p)];

wordSignals = [];
for i = 1:length(wordsplits)-2
    
    t = p(wordsplits(i):wordsplits(i+2)-1,:);
    wordSignals = cat(3, wordSignals, t);
    
end

clear p;

overSubjects = squeeze(mean(wordSignals,2));
overWords = squeeze(mean(wordSignals,3));

load('../../stimuli/mainVidInfo.mat')

wordInfo.words     = mainVidInfo.wordsShown;
wordInfo.present   = logical(mainVidInfo.timeSaid);
wordInfo.timeSaid = mainVidInfo.timeSaid;
wordInfo.timeShown = [5:5:99*5];
wordInfo.easyWords = [];
wordInfo.hardWords = [];

splitSeconds = [0:(pupilLength/fsEye/8):pupilLength/fsEye];

% split by time shown

for j = 1:length(splitSeconds)-1
    
    if(mod(j,2))
        t = find(wordInfo.timeShown > splitSeconds(j) & wordInfo.timeShown < splitSeconds(j+1));
        wordInfo.easyWords = [wordInfo.easyWords t];
    else
        t = find(wordInfo.timeShown > splitSeconds(j) & wordInfo.timeShown < splitSeconds(j+1));
        wordInfo.hardWords = [wordInfo.hardWords t];
    end
    
end

easySignals = wordSignals(:,:,wordInfo.easyWords);
hardSignals = wordSignals(:,:,wordInfo.hardWords);

% split by time said

if(timeSaid)

    fiveSecs = 5*fsEye;

    easySignals = [];
    hardSignals = [];

    for j = 1:length(wordInfo.timeSaid)

        if(wordInfo.timeSaid(j) == 0)
            continue;
        end

        section = [(wordInfo.timeSaid(j)*fsEye)-fiveSecs:(wordInfo.timeSaid(j)*fsEye)+fiveSecs-1];

        if(find(wordInfo.easyWords == j))

            easySignals = cat(3, easySignals, ppct(section,:));

        else

            hardSignals = cat(3, hardSignals, ppct(section,:));

        end

    end

end

% split by time said AND HIT OR MISS

if(timeSaid & separateHits)

    fiveSecs = 5*fsEye;

    for i = 1:length(subjectData)

        subjectData(i).easySignals = [];
        subjectData(i).hardSignals = [];
        subjectData(i).missedEasySignals = [];
        subjectData(i).missedHardSignals = [];

        for j = 1:length(wordInfo.timeSaid)

            if(wordInfo.timeSaid(j) == 0)
                continue;
            end

            section = [(wordInfo.timeSaid(j)*fsEye)-fiveSecs:(wordInfo.timeSaid(j)*fsEye)+fiveSecs-1];

            if(subjectData(i).missedWords(j))

                if(find(wordInfo.easyWords == j))

                    subjectData(i).missedEasySignals(:,end+1) = ppct(section,i);
    
                else
    
                    subjectData(i).missedHardSignals(:,end+1) = ppct(section,i);
    
                end

                continue;
            end

            if(find(wordInfo.easyWords == j))

                subjectData(i).easySignals(:,end+1) = ppct(section,i);

            else

                subjectData(i).hardSignals(:,end+1) = ppct(section,i);

            end

        end

        subjectData(i).easySignalsMean = mean(subjectData(i).easySignals,2);

        subjectData(i).hardSignalsMean = mean(subjectData(i).hardSignals,2);

        subjectData(i).missedEasySignalsMean = mean(subjectData(i).missedEasySignals,2);

        subjectData(i).missedHardSignalsMean = mean(subjectData(i).missedHardSignals,2);

    end

end

% easy average by subject HITS ONLY

figure;

e = cat(2,subjectData.easySignalsMean);

ea = mean(e, [1]);
ea = repmat(ea,[5000 1]);

plot([1:length(wordSignals)]/500, e-ea); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all HIT easy words for each subject');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeSaid_hitOnly_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeShown_hitOnly_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% easy average by subject MISSES ONLY

figure;

e = cat(2,subjectData.missedEasySignalsMean);

ea = mean(e, [1]);
ea = repmat(ea,[5000 1]);

plot([1:length(wordSignals)]/500, e-ea); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all MISS easy words for each subject');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeSaid_missOnly_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeShown_missOnly_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% easy average by subject HITS MINUS MISSES

figure;

e = cat(2,subjectData.missedEasySignalsMean);

ea = mean(e, [1]);
ea = repmat(ea,[5000 1]);

e = e-ea;

e2 = cat(2,subjectData.easySignalsMean);

ea2 = mean(e2, [1]);
ea2 = repmat(ea2,[5000 1]);

e2 = e2-ea2;

plot([1:length(wordSignals)]/500, e2-e); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all easy words (hits - misses)');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeSaid_hitMinusMiss_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeShown_hitMinusMiss_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% easy average by subject

figure;

plot([1:length(wordSignals)]/500, mean(easySignals,3)); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all easy words for each subject');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/easyWordAverageBySubject_timeShown_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% hard average by subject

figure;

plot([1:length(wordSignals)]/500, mean(hardSignals,3)); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all hard words for each subject');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/hardWordAverageBySubject_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/hardWordAverageBySubject_timeShown_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% easy grand average

figure;

shadedErrorBar([1:length(wordSignals)]/500, mean(easySignals,[2 3]), std(easySignals, 1, [2 3])/sqrt(length(subjects))); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all easy words across subjects (1/2 std err)');

set(gcf, 'position', [2146 160 797 707]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/easyWordGrandAverage_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/easyWordGrandAverage_timeShown_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% hard grand average

figure;

shadedErrorBar([1:length(wordSignals)]/500, mean(hardSignals,[2 3]), std(hardSignals, 1, [2 3])/sqrt(length(subjects))); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all hard words across subjects (std err)');

set(gcf, 'position', [2146 160 797 707]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/hardWordGrandAverage_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/hardWordGrandAverage_timeShown_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

close all;

% plot both together

figure;

e=shadedErrorBar([1:length(wordSignals)]/500, mean(easySignals,[2 3]), std(easySignals, 1, [2 3])/sqrt(length(subjects)),'lineProps',{'color',colors(2,:)}); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all words across subjects (std err)');

set(gcf, 'position', [2371 313 679 577]);

h=shadedErrorBar([1:length(wordSignals)]/500, mean(hardSignals,[2 3]), std(hardSignals, 1, [2 3])/sqrt(length(subjects)),'lineProps',{'color', colors(1,:)});

legend([e.mainLine h.mainLine], {'easy','hard'});

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/allWordGrandAverage_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/allWordGrandAverage_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

close all;

% subtract means from subject plots

% easy average by subject

figure;

ea = mean(easySignals, [1 3]);
ea = repmat(ea,[5000 1]);

easystd = std(mean(easySignals,3),1,2);

plot([1:length(wordSignals)]/500, mean(easySignals,3)-ea); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all easy words for each subject (minus mean)');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/easyWordAverageBySubjectMinusMean_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/easyWordAverageBySubjectMinusMean_timeShown_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

% hard average by subject

figure;

ha = mean(hardSignals, [1 3]);
ha = repmat(ha,[5000 1]);

plot([1:length(wordSignals)]/500, mean(hardSignals,3)-ha); hold on;
xline(5);

xticks([0 5 10]);
xticklabels({'-5','0','+5'});

xlabel('Seconds (s)');
ylabel('Pupil Dilation (% of baseline)');
title('Average of all hard words for each subject (minus mean)');

set(gcf, 'position', [2371 313 679 577]);

if(saveFlag)
    if(timeSaid)
        saveas(gcf, ['../../images/hardWordAverageBySubjectMinusMean_timeSaid_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    else
        saveas(gcf, ['../../images/hardWordAverageBySubjectMinusMean_timeShown_' char(datetime('today', 'format', 'MMddyy')) '.png']);
    end
end

close all;

% comparison with behavioral data

load('behavioralData.mat', 't');

toRemove = setdiff([t.subjectNumber], [subjects.subjectNumber]);

% remove non-matching subjects

for i = toRemove

    t(find([t.subjectNumber] == i)) = [];

end

% add behavior data to subject structure

for i = 1:length(t)

    subjects(i).f1 = t(i).f1;

end

e = mean(easyMean,2);
h = mean(hardMean,2);

figure;

for i = 1:length(subjects)

    p1 = plot(subjects(i).f1(1), e(i),'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;

    p2 = plot(subjects(i).f1(2), h(i),'color', colors(i,:), 'marker', 'x','linewidth',mywidth); hold on;

end

xlim([0.5 4.5]);
xticks([1.5 3.5]);
xticklabels({'-3dB', '-9dB'});
xlabel('Condition');

