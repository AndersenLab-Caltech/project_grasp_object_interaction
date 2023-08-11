%% script that checks whether the NSPs are swapped
% it compares 20210601 and 20210603 against 20210528 (the session that we know
% that NSPs are not swapped), to see whether some NSPs are swapped.
% It generates scatter plots comparing the median firing rate of each
% channel between two sessions - the points should be distributed along y=x
% if the NSPs are not swapped
​
%%
subj = hst.Subject('S2');
threshold_level = -4.5;
dateList = {'20210628', '20210603', '20210601'};
dateOfInterest = 1;
arrayOfInterest = [1,2];
arrayList = {'SMG_AIP', 'PMV', 'S1X_S1'};
% session_tue = hst.Session('20210601', subj);
% session_thu = hst.Session('20210603', subj);
% session_0528 = hst.Session('20210528', subj);
​
sessionList = cellfun(@(x) hst.Session(x, subj), dateList, 'uni', false);
​
%%
unitSummaryList = cellfun(@(x) x.unitSummary('unsorted', 'bad', 'FILE_SEARCH_ARGS', {'threshold', threshold_level}), sessionList, 'uni', false);
%%
FRList = cell(3,3);   % sessions * arrays
​
%%
for sessionIdx = 1:length(dateList)
    for i = 1:3
        FRList{sessionIdx, i} = zeros(1, 96);
        
        array = arrayList{i};
        
        for ch = 1:96
            temp = unitSummaryList{sessionIdx}{strcmpi(unitSummaryList{sessionIdx}{:,'nsp'}, array) & unitSummaryList{sessionIdx}{:,'channel'}==ch, 'firing_rate'};
            if ~isempty(temp), FRList{sessionIdx, i}(ch) = median(temp); end
        end
    end
end
​
​
%%
FRmax = ceil(max(cat(2, FRList{:})));
for arrayIdx = 1:length(arrayList)    
    for sessionIdx = 1:length(dateList)
        if sessionIdx == dateOfInterest, continue; end
        
        figure; scatter(FRList{dateOfInterest,arrayIdx}, FRList{sessionIdx,arrayIdx}, '.')
        title(sprintf('%s: %s vs. %s', arrayList{arrayIdx}, dateList{dateOfInterest}, dateList{sessionIdx}));  
        xlabel(sprintf('FR on %s', dateList{dateOfInterest})), ylabel(sprintf('FR on %s', dateList{sessionIdx}));
        xlim([0, FRmax ]), ylim([0, FRmax ]);
    end
end
% figure; scatter(FRList{3,1}, FRList{1,1}, '.'), title(sprintf('%s: %s vs. %s');  xlabel('FR on 0528'), ylabel('FR on 0601');     xlim([0, 90]), ylim([0, 90]);
% figure; scatter(FRList{3,1}, FRList{2,1}, '.'), title('SMG: 0528 vs. 0603');  xlabel('FR on 0528'), ylabel('FR on 0603');     xlim([0, 90]), ylim([0, 90]);
% figure; scatter(FRList{3,3}, FRList{2,3}, '.'), title('S1: 0528 vs. 0603');   xlabel('FR on 0528'), ylabel('FR on 0603');     xlim([0, 90]), ylim([0, 90]);
% figure; scatter(FRList{3,3}, FRList{1,3}, '.'), title('S1: 0528 vs. 0601');   xlabel('FR on 0528'), ylabel('FR on 0601');     xlim([0, 90]), ylim([0, 90]);
% 
% figure; scatter(FR_3{2}, FR_2{2}, '.'), title('PMv: 0528 vs. 0603');   xlabel('FR on 0528'), ylabel('FR on 0603');     xlim([0, 90]), ylim([0, 90]);
% figure; scatter(FR_3{2}, FRList{2}, '.'), title('PMv: 0528 vs. 0601');   xlabel('FR on 0528'), ylabel('FR on 0601');     xlim([0, 90]), ylim([0, 90]);