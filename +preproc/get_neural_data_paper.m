function [data_info] = get_neural_data_paper(task,varargin)


% Default variables values

%Control which brain regions to extract: 1: NSP1, 2: NSP2, 3:NSP3 and all combinations. Default SMG and PMV
[varargin,data_info.brain_region] = Utilities.argkeyval('brainregion',varargin, [1,2]); 
%Control which data is extracted: sorting (sorted spikes), noisy (includes noise units), unsorted (unsorted data), smoothed (halfkerner = 0 and causal = 0, how much the neuronal data is smoothed). 
[varargin,data_info.spikes_sorting] = Utilities.argkeyval('spikes',varargin, 'sorting');
%Control which trials are taken into account. Default: all trials
[varargin,data_info.trials] = Utilities.argkeyval('trials',varargin, 1:task.numTrials);
%Does not remove firing rates that are too low
[varargin,data_info.ratefilt_value] = Utilities.argkeyval('ratefilt',varargin,true); 
[varargin,min_timebin_length] = Utilities.argkeyval('min_timebin_length',varargin,158); 

Utilities.argempty(varargin);
%Save the name of the extracted arryas properly based on the subject

%Extract data depending on selected spike_sorting way
data_info.region = 'All';

data_info = preproc.spike_sorting_params(data_info);

disp(data_info.spikes)

debug = Debug.Debugger('testing');

%fr = firing rate
%relt: recording time. If <0, recorded before trial started.
%[fr,relt,featdef] = proc.task.bin(task,data_info.params,debug);

data_info.params.dt.cacheread = 0; %put cacheread to 0 to not automatically extract the S1 data again. 
data_info.params.tm.bufferpost = 1;
data_info.params.tm.bufferpre = 1;
% extract firing rates for trials. 
%[fr1,relt1,featdef1] = proc.task.bin(task,data_info.params,debug, 'trialidx', data_info.trials);
data_info.params.tm.bufferpost = 0; % do not extract data after trial.
data_info.params.tm.bufferpre = 0; % do not extract data before trial.
[fr,relt,featdef] = proc.task.bin(task,data_info.params,debug, 'trialidx', data_info.trials);

%NOTE: depending on how much data we take to extract, featdef can have
%different number of neurons because of the noratefilt kicking out process.
%E.g. I think that when we have a 1s buffer before and after, it kicks out
%more because we are in the ITI phase. 
%This explains why I am not extracting the same numeber of neurons in "old
%method" vs. "new method"

%realized issue: sometimes 1 trial is too short -> process it OUTSIDE so
%that I can pat that 1 trial with zeros instead of making all trials
%shorter. May or may not be in this dataset, so can ignore it 


%calculate phase labels 
phaseTimes = task.phaseTimes; % have to remove additional trials later otherwise calculations do not work out

%Calculate the mean time of each phase based on all the trials.
phaseTimesAvg = mean(phaseTimes - phaseTimes(:,1), 1);

numPhases = task.numPhases;
numTrials = task.numTrials; 
phaseTimeDuration = nan*zeros(numTrials,numPhases);

%calculate the duration of each phase
for n_phase = 1:(numPhases -1)
	phaseTimeDuration(:,n_phase) = -(phaseTimes(:,n_phase) -phaseTimes(:,n_phase+1));
end 

%calculate the duration of the last phase
for n_trials = 1:(numTrials-1)
    phaseTimeDuration(n_trials,numPhases) = -(phaseTimes(n_trials,numPhases) - phaseTimes(n_trials +1,1) );
end 
phaseTimeDuration(numTrials,numPhases) = nanmedian(phaseTimeDuration(:,numPhases));
medianPhaseTimeDuration = median(phaseTimeDuration);

if nnz(std(phaseTimeDuration) > data_info.params.spk.binwidth)
    disp('Big time difference in trial duration'); 
    %keyboard; 
end


%average phase duraction
avgPhaseTimeDuration = mean(phaseTimeDuration);
minPhaseTimeDuraction = min(phaseTimeDuration); 
%take shorted Action phase trial to not include ITI data
%phaseTimesAvg(end+1) = phaseTimesAvg(end) + avgPhaseTimeDuration(end);
phaseTimesAvg(end+1) = 7.9;

%calculate phase labels
phase_labels = arrayfun(@(x)phaseTimesAvg(x)<=relt, 1:length(phaseTimesAvg), 'UniformOutput', false);
phase_labels =  sum(cell2mat(phase_labels),2)';

try %test if 
    phase_times_idx(end+1) = find(relt > (medianPhaseTimeDuration(end) + phaseTimesAvg(end-1)),1);
    disp('Action phase trial duration is too short in one trial'); 
    %keyboard
catch 
    %minimum length of trial
    warning('Action phase shorter than expected: calculate relt and fr again'); 
    
    % pat shorter trials with nans instead of shortening all trials; 
    [fr,relt,featdef] = proc.task.bin(task,data_info.params,debug, 'UniformOutput', false,'trialidx', data_info.trials);
    
    %calculate size of trial:
    min_time_relt = cell2mat(cellfun(@(x) size(x,1), fr, 'UniformOutput', false));
    %find trials that are too short: 
    short_trials = find(min_time_relt < min_timebin_length);
    %find trials that are not too short: 
    ok_trials =setdiff(1:size(fr,1), short_trials);
    
    %pat the short trials with zeros
    for kk = 1:length(short_trials)
        s_tmp= short_trials(kk);
        length_dif = min_timebin_length - min_time_relt(s_tmp);
        %lengthen shorter trials
        fr{s_tmp} = [fr{s_tmp}; nan*zeros(length_dif, size(fr{s_tmp},2))];       
    end 
        
    %uniform relt and fr
    fr = cellfun(@(x) x(1:min_timebin_length,:), fr, 'UniformOutput', false);
    fr_tmp = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),fr,'UniformOutput',false));
    fr = permute(fr_tmp, [2,3, 1]);
    
    relt = relt{ok_trials(1)}(1:min_timebin_length);
    
    %re calculate phase labels
    phase_labels = arrayfun(@(x)phaseTimesAvg(x)<=relt, 1:length(phaseTimesAvg), 'UniformOutput', false);
    phase_labels =  sum(cell2mat(phase_labels),2)';

   % keyboard
    disp([' Number of trials needing to be longer: ' num2str(length(short_trials))]);
end 

%time bins to include in the analysis 
timebin_idx = ismember(phase_labels, 1:numPhases);

if size(timebin_idx,2) > min_timebin_length
    disp('Dataset has too many timebins - investigate')
    keyboard
end 

% extract relevant timebins and features
phase_labels = phase_labels(timebin_idx);
relt = relt(timebin_idx);
features_i = 1:length(featdef.nsp); %extract all features
fr_adapted = fr(timebin_idx, features_i, :);

% calculate the indexes where we start new phase
phase_times_idx(1) = 1;
phase_times_idx(1,2:numPhases) = find(diff(phase_labels))+1;

%save info
data_info.extracted_features = features_i; %save extracted features
data_info.phase_time_idx = phase_times_idx; %save time idx
data_info.featdef = featdef(features_i,:); %save feature defintion for extracted features
data_info.relt = relt; 
data_info.avg_phase_time_duration = avgPhaseTimeDuration;
data_info.fr_adapted = fr_adapted;
data_info.phase_labels = phase_labels; 
data_info.numPhases = numPhases;
data_info.numTrials = length(data_info.trials); 



end

