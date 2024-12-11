function step1_mergeEMG(subj, groupName, side)

%% initialize
clearvars -except subj groupName side
close all; clc;
fs = string(filesep) + string(filesep);
fPath = string(pwd)+fs;

if ~exist('subj','var') || isempty(subj), subj = "PS16"; else, subj = string(subj); end
if ~exist('groupName','var') || isempty(groupName), groupName = "young"; else, groupName = string(groupName); end
if ~exist('side','var') || isempty(side), side = "left"; else, groupName = string(groupName); end % options: both, left, right

projCode = "PS";
mergedSetName = "allsteps";
prepID = "incr0";
trialTypes = ["LEI", "LME", "REI", "RME"];
EMGfreq = 1e3;
resampfreq = 128;
%% construct necessary paths and files & adding paths
if ispc
    p2l.root = "Z:\BRaIN\"; p2l.git = "C:\~git\";
elseif isunix
    p2l.root = "~/BRaIN/"; p2l.git = "~/~git/";
end
p2l.eeglab = p2l.git + fs + "eeglab2019_0" + fs;
p2l.eegRepo = p2l.root + "eeg" + fs + projCode + fs + "EEG" + fs;
p2l.EEGsets = p2l.eegRepo + subj + fs + "EEG_sets" + fs; % Where you want to save your .set files
p2l.emg = p2l.eegRepo + subj + fs + "EMG" + fs;
p2l.ICA = p2l.eegRepo + subj + fs + "ICA" + fs; % ICA folder containing the prep folders
p2l.prep = p2l.ICA + prepID + fs; % increment 0 directory
p2l.studies = p2l.eegRepo + "STUDIES" + fs;
addpath(genpath(fPath))
addpath(genpath(fPath+fs+"funcs"))
addpath(p2l.eeglab)
if ~exist("pop_multifit.m","file"), eeglab; close; end
rmpath(p2l.eeglab + "plugins\MPT\dependency\propertyGrid\") % contains a faulty strjoin.m that crashes MATLAB

for t = trialTypes % create file paths to each trial's BDF file
    f2l.sets.(t) = subj + "_" + t + ".set";
    f2l.emg.(t) = subj + "_" + t + "_raw.mat";
end
f2l.comps = p2l.prep + subj + "_" + mergedSetName + "_ICA_STRUCT_rejbadchannels_diverse_select_comps.mat";
f2l.singlecomp = p2l.studies + "~clusters" + fs + groupName+"singleComponents.mat";

%% load the files
load(f2l.comps,"ICA_STRUCT")
singleComponents = load(f2l.singlecomp);
for t = trialTypes
   EEG =  pop_loadset( 'filename', char(f2l.sets.(t)), 'filepath', char(p2l.EEGsets));
   cprintf('Yellow',"loaded " + f2l.sets.(t) + "\n");

   [~, e63] = get_EEG_event_array(EEG,'63999','urevent'); % stepper events (square waves)
   [~, e64] = get_EEG_event_array(EEG,'64511','urevent'); % stepper events
   [~, e650] = get_EEG_event_array(EEG,'65023','urevent'); % stepper events (square waves)
   [~, e655] = get_EEG_event_array(EEG,'65535','urevent'); % stepper events
   EEG = pop_editeventvals(EEG,'delete',sort([e63 e64 e650 e655]));
   EEG=eeg_checkset(EEG,'eventconsistency', 'makeur');
   EEG=update_EEG(EEG,ICA_STRUCT,1,[],1,1);
   EEG.data = EEG.icaact(singleComponents.(subj),:);
   
   try
       EMG = loadEMG(p2l.emg+f2l.emg.(t),side);
   catch
       disp(subj+"_"+t+" EMG is not available")
   end
   % EMG highpass and rectify
   % We rectify because the negative singal does not have a biological
   % meaning.
   EMG.data = transpose(abs(highpass(EMG.data',1,EMGfreq)));
   % add the EEG first trigger time to EMG, so they have the same start
   EMG.time = EMG.time + EEG.trigger(2).latency/EEG.srate;
   % convert to milliseconds
   EMG.time = EMG.time * 1000;
   % resample EMG based on EEG time
   EMG.resampledData = transpose(interp1(EMG.time, EMG.data',EEG.times)); 
%    EMG.resampledData = EMG.resampledData * 1e6; % convert to micorvolt nit needed, let's do this using eeg_checkset 
   % add EMG to EEG
   EEG.data = [EEG.data;EMG.resampledData];
   EEG.muscleNames = EMG.muscleNames;
   EEG.icaweights = eye(size(EEG.data,1));
   EEG.icasphere = EEG.icaweights;
   EEG.icaact=[]; EEG.icawinv=[]; EEG.icachansind=[]; EEG.chanlocs=[];
   EEG.nbchan = size(EEG.data,1);
   EEG = eeg_checkset(EEG,'ica');
   % resample to a lower frequency
   EEG = pop_resample(EEG, resampfreq, 0.8, 0.4);
   EEG = eeg_checkset(EEG);
   
   % save raw EEG-EMG data
   EEG = pop_editset(EEG, 'setname', char(subj+"_"+t+"_ICEMG"));
   pop_saveset(EEG, 'filename', EEG.setname, 'filepath', char(p2l.EEGsets));
end


%% auxiliary functions
function EMG = loadEMG(EMGfile,side)
EMGstruct = load(EMGfile);
EMG = struct;
muscleNames = transpose(string(fieldnames(EMGstruct)));
EMG.muscleNames = [];
if side == "both"
    EMG.muscleNames = muscleNames;
elseif side == "left"
    for i = muscleNames, muscle = char(i); if strcmp('L',muscle(1)), EMG.muscleNames = [EMG.muscleNames i]; end, end
    EMG.muscleNames = [EMG.muscleNames "time"];
elseif side == "right"
    for i = muscleNames, muscle = char(i); if strcmp('R',muscle(1)), EMG.muscleNames = [EMG.muscleNames i]; end, end
    EMG.muscleNames = [EMG.muscleNames "time"];
end
for i = EMG.muscleNames
    if i == "time"
        EMG.time = reshape(EMGstruct.(i),1,[]);
    else
        EMG.data(i==EMG.muscleNames,:) = reshape(EMGstruct.(i),1,[]);
    end
end


