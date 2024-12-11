function step4_surrogate_analysis(subj, earlyFraction)
clearvars -except subj groupName earlyFraction;
close all; clc;
fs = string(filesep) + string(filesep);
fPath = string(pwd)+fs;
global epoch_timewindow
GUI_MODE            = 'nogui';  % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL     = 1;        % Verbosity Level (0=no/minimal output, 2=graphical output)
NumSample           = 200;      % number of samples for phase randomization                          
if ~exist('subj','var') || isempty(subj), subj = "PS27"; else, subj = string(subj); end
if ~exist('earlyFraction','var') || isempty(earlyFraction), earlyFraction = 1; end % fraction of pertubation to be flagged as early (or late)
if earlyFraction > .5, warning("early fraction is larger than half of the data, only one fraction (i.e., early) will be produced"); end
projCode = "PS";
trialTypes = ["LEI", "LME", "REI", "RME"];
if earlyFraction < .5
    phases = ["base_epoch_UPT" "early_epoch_CPT" "early_epoch_PPT" "late_epoch_CPT" ...
        "late_epoch_PPT" "post_epoch_UPT"];
else
    phases = ["early_epoch_CPT" "early_epoch_PPT"];
end
epoch_timewindow = [-1 1];
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
    p2l.surrogate.(t) = p2l.epoch.(t)+"surrogate"+fs;
    if ~isfolder(p2l.epoch.(t)), mkdir(p2l.epoch.(t)); end
    if ~isfolder(p2l.surrogate.(t)), mkdir(p2l.surrogate.(t)); end
    f2l.sets.(t) = subj+"_"+t+"_ICEMG" + ".set";
    for p = phases
        f2l.(p).(t).eeg = subj + "_" + t + "_" + p + "_CONN.set";
        f2l.(p).(t).dDTF = subj + "_" + t + "_" + p +'_ddtf08PConn.mat';
        f2l.(p).(t).sPConn = subj + "_" + t + "_" + p +'_SPConn.mat';
    end
end
%% compute the surrogate data
for t = trialTypes
    for p = phases
        try
            EEG = pop_loadset('filename', char(f2l.(p).(t).eeg), 'filepath', char(p2l.epoch.(t)));
            cprintf('Yellow',"loaded " + f2l.(p).(t).eeg + "\n");
            EEG.CAT.configs.stat_surrogateGen=struct([]);

            EEG = pop_stat_surrogateGen(EEG,GUI_MODE, ...
                'modelingApproach', EEG.CAT.configs.est_fitMVAR, ...
                'connectivityModeling',EEG.CAT.configs.est_mvarConnectivity, ...
                'mode',{'PhaseRand' 'nperms' NumSample}, ...
                'autosave',[], ...
                'verb',VERBOSITY_LEVEL);
            
            ddtf08PConn=EEG.CAT.PConn.dDTF08;
            save(p2l.surrogate.(t)+f2l.(p).(t).dDTF,'ddtf08PConn');
            SPConn=EEG.CAT.PConn.S;
            save(p2l.surrogate.(t)+f2l.(p).(t).sPConn,'SPConn');
            
            EEG = pop_editset(EEG, 'setname', [EEG.setname '_RANDPHASE']);
            pop_saveset(EEG, 'filename', EEG.setname, 'filepath', char(p2l.epoch.(t)));
        catch
            cprintf('Red',"There is a problem with " + t + "_" + p + "\n");
        end        
    end
end