% Do analysis of behavioral data for reward experiment 05/22
% Last edits: 05/16

clear all; close all; clc;

% windows paths
addpath(genpath('.'));
addpath('.');

% linux paths
addpath(genpath('.'));
addpath('.');

displayAnova = 'off';

maxDelay = 1.5;

behaviorLocation = '../../results/behavior/';
mywidth = 4; % width of lines to be drawn

numBlocks = 8;

t = table2struct(readtable(['../../subjects.xlsx']));

t = t(logical([t.include]));

load('.','splitVidInfo');

for i = 1:length(t)
    
    load([behaviorLocation t(i).behaviorFile '.mat'])
    
    t(i).behaviorResults = experimentData;

    t(i).conditionOrder = t(i).behaviorResults.conditionOrder;

    t(i).adjustedPresses = cellfun(@(x,kpt) x-kpt,t(i).behaviorResults.keyPressTimes, ...
        mat2cell(t(i).behaviorResults.soundOnsets',ones(1,8))', 'uniformoutput', false);

    for j = 1:length(t(i).behaviorResults.keyPressTimes)
    
        k = t(i).behaviorResults.keyPressTimes{j};

        t(i).timeSaid = splitVidInfo.timeSaid;

    end
    
    t(i).numPossible  = [0 0 0 0];
    t(i).numResponses = [0 0 0 0];
    t(i).numHits      = [0 0 0 0];
    t(i).numMisses    = [0 0 0 0];
    t(i).numFA        = [0 0 0 0];

    t(i).responseHit = cell(1,numBlocks);

    % find number possible for each subject/condition
    for j = 1:numBlocks

        t(i).numPossible(t(i).conditionOrder(j)) = t(i).numPossible(t(i).conditionOrder(j)) + ...
            length(find(t(i).timeSaid(:,j)));

        t(i).numResponses(t(i).conditionOrder(j)) = t(i).numResponses(t(i).conditionOrder(j)) + ...
            length(t(i).behaviorResults.keyPressTimes{j});

        t(i).responseHit{j} = zeros(1,length(t(i).behaviorResults.keyPressTimes{j}));

    end
    
    for j = 1:numBlocks % label each key press as hit or false alarm

        ts = t(i).timeSaid(:,j)';

        kpt = t(i).adjustedPresses{j};

        rh = t(i).responseHit{j};

        rt = t(i).responseHit{j};

        for k = 1:length(kpt)

            if(find(kpt(k) > ts & kpt(k) < (ts+maxDelay)))
                rh(k) = 1;
                rt(k) = min(abs(ts - kpt(k)));
            end

        end

        t(i).responseHit{j}  = rh;
        t(i).reactionTime{j} = rt;

    end

    h = cellfun(@(x) length(find(x)), t(i).responseHit);

    for j = 1:numBlocks % sum up hits

        t(i).numHits(t(i).conditionOrder(j)) = t(i).numHits(t(i).conditionOrder(j)) + h(j);

    end

    t(i).numMisses = t(i).numPossible  - t(i).numHits;
    t(i).numFA     = t(i).numResponses - t(i).numHits;
    
    t(i).f1 = t(i).numHits./(t(i).numHits + 0.5*(t(i).numFA + t(i).numMisses));

    % do f1 by block

    for j = 1:numBlocks

        t(i).numPossibleByBlock(j) = length(find(t(i).timeSaid(:,j))); % find non-zero time said to get # possible
        t(i).numHitsByBlock(j)     = length(find(t(i).responseHit{j}));
        t(i).numFAByBlock(j)       = length(t(i).responseHit{j}) - t(i).numHitsByBlock(j);
        t(i).numMissesByBlock(j)   = t(i).numPossibleByBlock(j) - t(i).numHitsByBlock(j);

    end

    t(i).f1ByBlock = t(i).numHitsByBlock./(t(i).numHitsByBlock + 0.5*(t(i).numFAByBlock + t(i).numMissesByBlock));
    
end

% ANOVA analysis

f1 = cat(1,t.f1);
subject = repmat([1:length(t)]', [1 4]);
reward  = repmat([1 1 0 0], [length(t) 1]);
noise   = repmat([3 9 3 9], [length(t) 1]);

behaviorTable = table(f1(:),subject(:), reward(:), noise(:), 'VariableNames',{'f1', 'subject', 'reward', 'noise'});

% lme = fitlme(behaviorTable, '')

[p tab stats terms] = anovan(f1(:), {subject(:), reward(:), noise(:)}, 'model', 'interaction', 'random', [1], ...
    'varnames', {'subject', 'reward', 'noise'}, 'display', displayAnova);

% plot behavior fig - reward connected

% 1 - 3db reward
% 2 - 9db reward
% 3 - 3db no
% 4 - 9db no

colors = cat(1, brewermap(NaN,'Set1'), brewermap(NaN,'Set2'), brewermap(NaN,'Set3'), brewermap(NaN,'Dark2'), brewermap(NaN,'Paired'), brewermap(NaN,'Accent'));
colors = cat(2,colors,repmat([0.7],length(colors),1));
colors = repmat(colors, 2, 1);

f = figure;

for i = 1:length(t)

    plot(1, t(i).f1(1), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(1.5, t(i).f1(3), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([1 1.5], t(i).f1([1 3]), 'color', colors(i,:),'linewidth',mywidth);

    plot(2, t(i).f1(2), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(2.5, t(i).f1(4), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([2 2.5], t(i).f1([2 4]), 'color', colors(i,:),'linewidth',mywidth);

end

means = mean(cat(1,t.f1));
errors = std(cat(1,t.f1))/sqrt(length(t));

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
xticklabels({'-3dB R', '-3dB N', '-9dB R', '-9dB N'})
ylabel('F1 Score')
title('Behavior Results')

f1 = cat(1,t.f1);

[h p] = ttest(f1(:,1),f1(:,3));

sigstar({[1.25 2.25] [2 2.5]}, [0.001, 0.05]);

set(gcf, 'position', [950 497 403 375]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);
box off;

saveas(gcf, ['../../images/behavioralPerformance_rewardConnected_' char(datetime('today', 'format', 'MMddyy')) '.png'])
saveas(gcf, ['../../images/behavioralPerformance_rewardConnected_' char(datetime('today', 'format', 'MMddyy')) '.epsc'])
saveas(gcf, ['../../images/behavioralPerformance_rewardConnected_' char(datetime('today', 'format', 'MMddyy')) '.fig'])

close;

% behavior figure - noise connected

f = figure;

for i = 1:length(t)

    plot(1, t(i).f1(1), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(1.5, t(i).f1(2), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([1 1.5], t(i).f1([1 2]), 'color', colors(i,:),'linewidth',mywidth);

    plot(2, t(i).f1(3), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth); hold on;
    plot(2.5, t(i).f1(4), 'color', colors(i,:), 'marker', 'o','linewidth',mywidth);
    plot([2 2.5], t(i).f1([3 4]), 'color', colors(i,:),'linewidth',mywidth);

end

xlim([0.7 2.8]);
xticks([1 1.5 2 2.5])
xticklabels({'-3dB R', '-9dB R', '-3dB N', '-9dB N'})
ylabel('F1 Score')
title('Behavioral performance')

f1 = cat(1,t.f1);

[h p] = ttest(f1(:,1),f1(:,2));

sigstar({[1 1.5]}, p)

[h p] = ttest(f1(:,3),f1(:,4));

sigstar({[2 2.5]}, p)

% ylim([0 1.1])

set(gcf, 'position', [950 497 403 375]);
set(findall(gcf,'-property','FontSize'),'FontSize',16);
box off;

saveas(gcf, ['../../images/behavioralPerformance_noiseConnected_' char(datetime('today', 'format', 'MMddyy')) '.png'])
saveas(gcf, ['../../images/behavioralPerformance_noiseConnected_' char(datetime('today', 'format', 'MMddyy')) '.epsc'])
saveas(gcf, ['../../images/behavioralPerformance_noiseConnected_' char(datetime('today', 'format', 'MMddyy')) '.fig'])

close;

% do reaction time plot

f = figure;

for i = 1:length(t)

    rtr = cat(2, t(i).reactionTime{find(t(i).conditionOrder == 1 | t(i).conditionOrder == 2)});
    rtr = rtr(find(rtr));
    rtn = cat(2, t(i).reactionTime{find(t(i).conditionOrder == 4 | t(i).conditionOrder == 3)});
    rtn = rtn(find(rtn));

    t(i).rtr = rtr;
    t(i).rtn = rtn;

end

% do histogram
hr = histogram(cat(2,t.rtr),15); hold on; hn = histogram(cat(2,t.rtn), hr.NumBins);

% adjust
h3.FaceColor = [55,126,184]/255;
h9.FaceColor = [228,26,28]/255;

legend({'-3dB R', '-9dB R', '-3dB N', '-9dB N'});

xlabel('Reaction Time (s)');
ylabel('Number of Responses');

set(findall(gcf,'-property','FontSize'),'FontSize',16);

% set(gcf, 'position', [425 547 414 339]);
set(gcf, 'position', [2371 532 418 358]);
box off;

% save
saveas(gcf, ['../../images/reactionTime_' char(datetime('today', 'format', 'MMddyy')) '.png'])
saveas(gcf, ['../../images/reactionTime_' char(datetime('today', 'format', 'MMddyy')) '.epsc'])
saveas(gcf, ['../../images/reactionTime_' char(datetime('today', 'format', 'MMddyy')) '.fig'])

close;

% generate behavioral "signals"

load('.')

fsb = 60;

behaviorLen = 224.5; % in seconds

bTime   = linspace(0, behaviorLen, fsb*behaviorLen);
bSignal = zeros(fsb*behaviorLen, numBlocks);

for i = 1:length(t)

    behavior = t(i).behaviorResults;

    missedWords = zeros(size(t(i).timeSaid));

    keyTimes = t(i).adjustedPresses;

    for j = 1:numBlocks

        for k = 1:length(missedWords)

            if(t(i).timeSaid(k,j) == 0)
                continue;
            end

            if(~any(find(t(i).timeSaid(k,j) < keyTimes{j} & t(i).timeSaid(k,j)+maxDelay > (keyTimes{j}))))
                missedWords(k,j) = 1;
            end

        end

    end

    t(i).missedWords = missedWords;
    t(i).bSignal = bSignal;
    t(i).bSignalPressed = bSignal;
    t(i).bSignalPressedHit = bSignal;
    t(i).bSignalPressedMiss = bSignal;

    for j = 1:length(t(i).missedWords)

        for k = 1:numBlocks

            if(t(i).timeSaid(j,k) == 0)
                continue;
            end

            [x, ix] = min(abs(bTime - t(i).timeSaid(j,k)));

            if(t(i).missedWords(j,k) == 0 & t(i).timeSaid(j,k) ~= 0)
                t(i).bSignal(ix,k) = 1;
            else
                t(i).bSignal(ix,k) = -1;
            end

        end

    end

    for j = 1:numBlocks

        for k = 1:length(t(i).adjustedPresses{j})

            [x, ix] = min(abs(bTime - t(i).adjustedPresses{j}(k)));

            t(i).bSignalPressed(ix, j) = 1;

        end

        tmp = t(i).adjustedPresses{j};
        tmp = tmp(logical(t(i).responseHit{j}));

        for k = 1:length(tmp)

            [x, ix] = min(abs(bTime - tmp(k)));

            t(i).bSignalPressedHit(ix, j) = 1;

        end

        tmp = t(i).adjustedPresses{j};
        tmp = tmp(~logical(t(i).responseHit{j}));

        for k = 1:length(tmp)

            [x, ix] = min(abs(bTime - tmp(k)));

            t(i).bSignalPressedMiss(ix, j) = 1;

        end

    end

end

save('missedWordData.mat', 't');

% find most missed words

totalMissedWords = {};
totalHitWords    = {};

for i = 1:length(t)

    totalMissedWords = [totalMissedWords; splitVidInfo.words(logical(t(i).missedWords))];

    totalHitWords = [totalHitWords; splitVidInfo.words(~logical(t(i).missedWords))];

end

[missedWordCounts] = countWords(totalMissedWords);

[totalWordCounts] = countWords(splitVidInfo.words(:));

% find words related to false alarms

for i = 1:length(t)

    t(i).falseAlarmWords = {};
    t(i).falseAlarmWords{1, end+1} = 'word';
    t(i).falseAlarmWords{2, end}   = 'section';
    t(i).falseAlarmWords{3, end}   = 'index';
    t(i).falseAlarmWords{4, end}   = 'time shown';
    t(i).falseAlarmWords{5, end}   = 'time press';

    for j = 1:8

        fat = t(i).adjustedPresses{j}(~t(i).responseHit{j});
        
        ti = [0:5:224.5];

        for k = 1:length(fat)

            [x, ix] = min(abs(ti - fat(k)));

            if(ti(ix) > fat(k))
                ix = ix - 1;
            end

            t(i).falseAlarmWords{1, end+1} = splitVidInfo.words{ix, j};
            t(i).falseAlarmWords{2, end}   = j;
            t(i).falseAlarmWords{3, end}   = ix;
            t(i).falseAlarmWords{4, end}   = ti(ix);
            t(i).falseAlarmWords{5, end}   = fat(k);

        end

    end

end

% count false alarm words

falseAlarmWords = {};

for i = 1:length(t)

    falseAlarmWords = [falseAlarmWords t(i).falseAlarmWords(1, 2:end)];

end