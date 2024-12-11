function step2_splitEpoch_preprocess(subj, earlyFraction)
clearvars -except subj groupName earlyFraction;
close all; clc;
fs = string(filesep) + string(filesep);
fPath = string(pwd)+fs;
global epoch_timewindow

if ~exist('subj','var') || isempty(subj), subj = "PS04"; else, subj = string(subj); end
if ~exist('earlyFraction','var') || isempty(earlyFraction), earlyFraction = 0.999; end % fraction of pertubation to be flagged as early (or late)
if earlyFraction > .5, warning("early fraction is larger than half of the data, only one fraction (i.e., early) will be produced"); end
projCode = "PS";
trialTypes = ["LEI", "LME", "REI", "RME"];
if earlyFraction > .5
    timeSect = ["base","early","post"]; %time sections
else
    timeSect = ["base","early","late","post"]; %time sections
end
epoch_timewindow = [-2 2];

%% construct necessary paths and files & adding paths
if ispc
    p2l.root = "Z:\BRaIN\"; p2l.git = "C:\~git\";
elseif isunix
    p2l.root = "~/BRaIN/"; p2l.git = "~/~git/";
end
p2l.eeglab = p2l.git + fs + "eeglab2019_0" + fs;
p2l.eegRepo = p2l.root + "eeg" + fs + projCode + fs + "EEG" + fs;
p2l.EEGsets = p2l.eegRepo + subj + fs + "EEG_sets" + fs; % Where you want to save your .set files

addpath(genpath(fPath))
addpath(genpath(fPath+fs+"funcs"))
addpath(p2l.eeglab)
if ~exist("pop_multifit.m","file"), eeglab; close; end
rmpath(p2l.eeglab + "plugins\MPT\dependency\propertyGrid\") % contains a faulty strjoin.m that crashes MATLAB

for t = trialTypes % create file paths to each trial's BDF file
    p2l.epoch.(t) = p2l.EEGsets + "epoch_fraction" +string(earlyFraction) + fs + "conn" + fs + t +fs;
    if ~isfolder(p2l.epoch.(t)), mkdir(p2l.epoch.(t)); end
    f2l.sets.(t) = subj+"_"+t+"_ICEMG" + ".set";
    for s = timeSect
        f2l.(s).(t) = subj + "_" + t + "_" + s + ".set";
    end
end

%% split
% here we don't need oultier recording because we want to analyze
% connectvity around the pertubations. Nevertheless, we need to do the
% outlier analysis to hve better estimate for perturbation latency for
% unperturbaed strdes
for t = trialTypes
    EEG = pop_loadset('filename', char(f2l.sets.(t)), 'filepath', char(p2l.EEGsets));
    cprintf('Yellow',"loaded " + f2l.sets.(t) + "\n");
    
    switch t
        case {"LEI", "LME"}
           [PPL.latencies, PPL.n] = get_EEG_event_array(EEG,'PPL','latency');
           [PPT.latencies, PPT.n] = get_EEG_event_array(EEG,'PPT','latency');
           [POR.latencies, POR.n] = get_EEG_event_array(EEG,'POR','latency');
           [UPL.latencies, UPL.n] = get_EEG_event_array(EEG,'UPL','latency');
           [UPT.latencies, UPT.n] = get_EEG_event_array(EEG,'UPT','latency');
           [CPL.latencies, CPL.n] = get_EEG_event_array(EEG,'CPL','latency');
           [CPT.latencies, CPT.n] = get_EEG_event_array(EEG,'CPT','latency');
    
           perturb_latency = [PPT.latencies{:}]-[PPL.latencies{:}];
           [latency_outlier,L,U,~] = isoutlier(perturb_latency,"mean","ThresholdFactor",3);
           ave_perturb_latency = round(mean(perturb_latency(~latency_outlier)));
           for i = 1:length(UPL.n)
                   EEG.event(UPT.n(i)).latency = UPL.latencies{i}+ave_perturb_latency;
           end
           for i = 1:length(CPL.n)
                   EEG.event(CPT.n(i)).latency = CPL.latencies{i}+ave_perturb_latency;
           end
           EEG=eeg_checkset(EEG,'eventconsistency', 'makeur');
           EEGbase = pop_select(EEG, 'point', [1 EEG.event(PPL.n(1)+5).latency]);
           EEGbase.setname = char(subj + "_" + t + "_base");
           EEGbase = eeg_checkset(EEGbase,'eventconsistency', 'makeur');
           pop_saveset(EEGbase,'filename', EEGbase.setname, 'filepath', char(p2l.epoch.(t)));
           
           PPL_early = ceil(length(PPL.n)*earlyFraction); % Changed to earlyFraction to able to control the fraction 10-07-19
           if earlyFraction < 0.5, PPL_late = length(PPL.n) - ceil(length(PPL.n)*earlyFraction); end
           EEGearly = pop_select(EEG, 'point', [EEG.event(PPL.n(1)-5).latency EEG.event(PPL.n(PPL_early)+5).latency]);
           EEGearly.setname = char(subj + "_" + t + "_early");
           EEGearly = eeg_checkset(EEGearly,'eventconsistency', 'makeur');
           pop_saveset(EEGearly,'filename', [EEGearly.setname], 'filepath', char(p2l.epoch.(t)));
           
           if earlyFraction < 0.5
               EEGlate = pop_select(EEG, 'point', [EEG.event(PPL.n(PPL_late)-5).latency EEG.event(POR.n(end)+5).latency]);
               EEGlate.setname = char(subj + "_" + t + "_late");
               EEGlate = eeg_checkset(EEGlate,'eventconsistency', 'makeur');
               pop_saveset(EEGlate,'filename', EEGlate.setname, 'filepath', char(p2l.epoch.(t)));
           end
           
           EEGpost = pop_select(EEG, 'point', [EEG.event(POR.n(end)-5).latency max([EEG.event.latency])]);
           EEGpost.setname = char(subj + "_" + t + "_post");
           EEGpost = eeg_checkset(EEGpost,'eventconsistency', 'makeur');
           pop_saveset(EEGpost,'filename', EEGpost.setname, 'filepath', char(p2l.epoch.(t)));

        case {"REI", "RME"}
           [PPR.latencies, PPR.n] = get_EEG_event_array(EEG,'PPR','latency');
           [PPT.latencies, PPT.n] = get_EEG_event_array(EEG,'PPT','latency');
           [POL.latencies, POL.n] = get_EEG_event_array(EEG,'POL','latency');
           [UPR.latencies, UPR.n] = get_EEG_event_array(EEG,'UPR','latency');
           [UPT.latencies, UPT.n] = get_EEG_event_array(EEG,'UPT','latency');
           [CPR.latencies, CPR.n] = get_EEG_event_array(EEG,'CPR','latency');
           [CPT.latencies, CPT.n] = get_EEG_event_array(EEG,'CPT','latency');

           perturb_latency = [PPT.latencies{:}]-[PPR.latencies{:}];
           [latency_outlier,L,U,~] = isoutlier(perturb_latency,"mean","ThresholdFactor",3);
           ave_perturb_latency = round(mean(perturb_latency(~latency_outlier)));
           for i = 1:length(UPR.n)
                   EEG.event(UPT.n(i)).latency = UPR.latencies{i}+ave_perturb_latency;
           end
           for i = 1:length(CPR.n)
                   EEG.event(CPT.n(i)).latency = CPR.latencies{i}+ave_perturb_latency;
           end
           EEG = eeg_checkset(EEG,'eventconsistency', 'makeur');

           EEGbase = pop_select(EEG, 'point', [1 EEG.event(PPR.n(1)+5).latency]);
           EEGbase.setname = char(subj + "_" + t + "_base");
           EEGbase = eeg_checkset(EEGbase,'eventconsistency', 'makeur');
           pop_saveset(EEGbase,'filename', EEGbase.setname, 'filepath', char(p2l.epoch.(t)));
           
           PPR_early = ceil(length(PPR.n)*earlyFraction); % Changed to earlyFraction to able to control the fraction 10-07-19
           if earlyFraction < 0.5, PPR_late = length(PPR.n) - ceil(length(PPR.n)*earlyFraction); end
           EEGearly = pop_select(EEG, 'point', [EEG.event(PPR.n(1)-5).latency EEG.event(PPR.n(PPR_early)+5).latency]);
           EEGearly.setname = char(subj + "_" + t + "_early");
           EEGearly = eeg_checkset(EEGearly,'eventconsistency', 'makeur');
           pop_saveset(EEGearly,'filename', EEGearly.setname, 'filepath', char(p2l.epoch.(t)));
           
           if earlyFraction < 0.5
               EEGlate = pop_select(EEG, 'point', [EEG.event(PPR.n(PPR_late)-5).latency EEG.event(POL.n(end)+5).latency]);
               EEGlate.setname = char(subj + "_" + t + "_late");
               EEGlate = eeg_checkset(EEGlate,'eventconsistency', 'makeur');
               pop_saveset(EEGlate,'filename', EEGlate.setname, 'filepath', char(p2l.epoch.(t)));
           end
           
           [EEGpost, eventidx] = pop_select(EEG, 'point', [EEG.event(POL.n(end)-5).latency max([EEG.event.latency])]);
           EEGpost.setname = char(subj + "_" + t + "_post");
           EEGpost = eeg_checkset(EEGpost,'eventconsistency', 'makeur');
           pop_saveset(EEGpost,'filename', EEGpost.setname, 'filepath', char(p2l.epoch.(t)));
    end
end

%% epoch
% This is a spacial case of run6, where there is ony one event to epoch to,
% plus/minus X millisecond. So, no warping is required. Basically, I need to
% create pre, post, perturb and catches with earlyFraction in mind.
for t = trialTypes
    for s = timeSect
        EEG = pop_loadset('filename', char(f2l.(s).(t)), 'filepath', char(p2l.epoch.(t)));
        cprintf('Yellow',"loaded " + f2l.(s).(t) + "\n");
        switch s
            case {"early", "late"}
                for i = ["PPT", "CPT"]
                    [EEGep.(t).(s).(i), ~] = create_epoch(EEG,i,[]);
                    [EEGep.(t).(s).(i), EEGep.(t).(s).(i).modelOrder] = preProcess_eeg(EEGep.(t).(s).(i),0);
                    pop_saveset(EEGep.(t).(s).(i), 'filename', EEGep.(t).(s).(i).setname, 'filepath', char(p2l.epoch.(t)));
                end
            case {"base", "post"}
                [EEGep.(t).(s).UPT, ~] = create_epoch(EEG,"UPT",[]);
                [EEGep.(t).(s).UPT, EEGep.(t).(s).UPT.modelOrder] = preProcess_eeg(EEGep.(t).(s).UPT,0);
                pop_saveset(EEGep.(t).(s).UPT, 'filename', EEGep.(t).(s).UPT.setname, 'filepath', char(p2l.epoch.(t)));
        end
    end
end

function [EEG,e_seq] = create_epoch(EEG,e_event,e_seq, randSelect)
global epoch_timewindow
if ~exist('randSelect','var') || isempty(randSelect), randSelect = 0; end
EEG = pop_epoch( EEG, {char(e_event)}, epoch_timewindow,...
    'newname', [EEG.setname '_epoch_' char(e_event)], 'epochinfo', 'yes');
EEG = eeg_checkset(EEG,'eventconsistency', 'makeur');
if randSelect
    EEG = pop_select(EEG,'trial',sort(randi(length(EEG.epoch),1,randSelect)));
    EEG = eeg_checkset(EEG,'eventconsistency', 'makeur');
end
e_seq = [e_event e_seq];

function [EEG,modelOrder] = preProcess_eeg(EEG,identifyMinOrder)
if ~exist('identifyMinOrder','var') || isempty(identifyMinOrder), identifyMinOrder = 0; end
modelOrder = 32; % default value, based on Steve Peterson's paper, to modify, change identifyMinOrder to 1
components          = [];                           % these are the components/channels to which we'll fit our multivariate model (if empty, select all in ica weight matrix)
windowLengthSec     = 1.3;                          % sliding window length in seconds
windowStepSizeSec   = 0.01;                         % sliding window step size in seconds
GUI_MODE            = 'nogui';                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL     = 0;                            % Verbosity Level (0=no/minimal output, 2=graphical output)
modelOrderCrit='elbow'; %'min';
if isempty(components), components=1:size(EEG.icaweights,1);end
% convert list of components to cell array of strings
componentNames = strtrim(cellstr(num2str(components')));

[EEG] = pop_pre_prepData(EEG,GUI_MODE, ...
    'VerbosityLevel',VERBOSITY_LEVEL,   ...
    'SignalType',{'Components'},  ...
    'VariableNames',componentNames,   ...
    'Detrend',  ...
    {'verb' VERBOSITY_LEVEL ...
    'method' {'linear'} ...
    'piecewise' ...
    {'seglength' 0.33   ...
    'stepsize' 0.0825} ...
    'plot' false},  ...
    'NormalizeData',    ...
    {'verb' 0       ...
    'method' {'time' 'ensemble'}},   ...
    'resetConfigs',true,    ...
    'badsegments',[],       ...
    'newtrials',[],         ...
    'equalizetrials',false);

if ~identifyMinOrder
    GUI_MODE = 'cfg_only';
end
[EEG,cfg] = pop_est_selModelOrder(EEG,GUI_MODE, ...
    'modelingApproach',         ...
    {'Segmentation VAR'     ...
    'algorithm' {'Vieira-Morf'} ... % Ridge Regression  Vieira-Morf(defualt)  'Group Lasso (ADMM)','GPU','on'
    'winStartIdx' []    ...
    'winlen'  windowLengthSec    ...
    'winstep' windowStepSizeSec  ...
    'taperfcn' 'blackmanharris'    ...
    'epochTimeLims' []      ...
    'prctWinToSample' 100   ...
    'normalize' [] ...
    'detrend' {'method' 'constant'} ...
    'verb' VERBOSITY_LEVEL},      ...
    'morderRange',[1 40], ... %30] ,  ...
    'downdate',true,        ...
    'runPll',[],            ...
    'icselector',{'sbc' 'aic' 'fpe' 'hq'},  ...
    'winStartIdx',[],       ...
    'epochTimeLims',[],     ...
    'prctWinToSample',80,   ...
    'plot', [], ...
    'verb',VERBOSITY_LEVEL);
if identifyMinOrder
     if strcmp(modelOrderCrit,'min')
         modelOrder = ceil(mean(EEG(1).CAT.IC.hq.popt));
     elseif strcmp(modelOrderCrit,'elbow')
         modelOrder = ceil(mean(EEG(1).CAT.IC.aic.pelbow));
     end
else
    EEG.CAT.configs.est_selModelOrder = cfg;
end
disp(['Model Order Selected is: ' num2str(modelOrder)]);

