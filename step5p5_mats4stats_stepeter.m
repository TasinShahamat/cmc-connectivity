function step5p5_mats4stats(groupName, earlyFraction)
clearvars -except subj groupName earlyFraction;
close all; clc;
fs = string(filesep) + string(filesep);
fPath = string(pwd)+fs;
global timewindow
timewindow = [-0.4,-0.1];
alpha = 0.01; % which one, 0.05 or 0.01 ?
windowLengthSec     = 0.4;                          % sliding window length in seconds
windowStepSizeSec   = 0.02;                         % sliding window step size in seconds
GUI_MODE            = 'nogui';                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL     = 0;                            % Verbosity Level (0=no/minimal output, 2=graphical output)
if ~exist('groupName','var') || isempty(groupName), groupName = "young"; else, groupName = string(groupName); end
if ~exist('earlyFraction','var') || isempty(earlyFraction), earlyFraction = 0.999; end % fraction of pertubation to be flagged as early (or late)
if earlyFraction > .5, warning("early fraction is larger than half of the data, only one fraction (i.e., early) will be produced"); end
if ~exist('subjs','var') || isempty(subjs)
    if groupName == "young"
        subjs = ["PS04" "PS05" "PS06" "PS07" "PS15" "PS16" "PS17" "PS18" "PS19" ...
            "PS20" "PS21" "PS22" "PS23" "PS24" "PS25" "PS26" "PS27"];
    elseif groupName == "old"
        subjs = ["PS08" "PS09" "PS10" "PS28" "PS29" "PS31" "PS32" "PS33" ...
            "PS34" "PS35" "PS36"];
    end
else, subjs = string(subjs); % 17 subjs
end
projCode = "PS";
% trialTypes = ["LEI", "LME", "REI", "RME"];
trialTypes = ["LEI", "LME"];  % left side muscles are considered. LME is showing strange stuff suddenly (only showing dDTF values along ppc_r). Need tp debug
if earlyFraction < .5
    phases = ["early_epoch_CPT" "early_epoch_PPT" "late_epoch_CPT" ...
        "late_epoch_PPT"];
else
%         phases = ["early_epoch_CPT" "early_epoch_PPT"]; % just running connectivty for perturbed and catch
    phases = ["early_epoch_PPT"];
end

% The muscle and cluster orders are critically important, becuase the ICs are sorted in this order
muscles = ["LTA","LSO","LRF","LST","LAD","LPD"]; % look in step1 EMG struct for the order ,"RTA","RSO","RRF","RST","RAD","RPD"
if groupName == "young"
    clustNames.cls04="ACC"; clustNames.cls11="SMA_r";clustNames.cls15="SMA_l";clustNames.cls08="PPC_r";clustNames.cls16="PPC_l";
    bNames = ["ACC", "SMA_r", "SMA_l", "PPC_r", "PPC_l"];
else
    clustNames.cls07="Motor"; clustNames.cls03="PPC_r";clustNames.cls06="PPC_l";
    bNames = ["Motor", "PPC_r", "PPC_l"];
end
numComp = length(muscles) + length(bNames);
%%
if ispc
    p2l.root = "Z:\BRaIN\"; p2l.git = "C:\~git\";
elseif isunix
    p2l.root = "~/BRaIN/"; p2l.git = "~/~git/";
end
p2l.eeglab = "D:\eeglab2019_0_seyed" + fs;
p2l.eegRepo = p2l.root + "eeg" + fs + projCode + fs + "EEG" + fs;
% p2l.EEGsets = p2l.eegRepo + subj + fs + "EEG_sets" + fs; % Where you want to save your .set files
p2l.studies = p2l.eegRepo + "STUDIES" + fs;
p2l.studyConn = p2l.studies + fs + "conn" + fs;
if ~exist(p2l.studyConn,"file"), mkdir(p2l.studyConn); end
addpath(genpath(fPath))
addpath(genpath(fPath+fs+"funcs"))
addpath(p2l.eeglab)
if ~exist("pop_multifit.m","file"), eeglab; close; end
rmpath(p2l.eeglab + "plugins\MPT\dependency\propertyGrid\") % contains a faulty strjoin.m that crashes MATLAB
for s = subjs
    for t = trialTypes
        p2l.conn.(s).(t) = p2l.eegRepo+s+fs+"EEG_sets"+fs+"epoch_fraction"+string(earlyFraction)+fs+"conn"+fs+t+fs;
        for p = phases
            f2l.conn.(s).(t).(p) = s+"_"+t+"_"+p+"_CONN.set"; % need only file name
        end
    end
end

f2l.singlecomp = p2l.studies + "~clusters" + fs + groupName+"singleComponents.mat";

%% construct the big connectivty array
singClust = load(f2l.singlecomp); % conains both a list of comps per subject, and list of clusters and their comps
for t = trialTypes
    for p = phases
        connStruct_boot = cell(numComp,numComp); % dDTF08 length(muscles)+length(fielnames(clusterNames)) is numComp
        connStructS_boot = connStruct_boot;
        numConns_boot = zeros(numComp,numComp); % this is to track which clusters had what comps.
        EEG = [];
        for s = subjs
            try
                EEG = pop_loadset('filename',char(f2l.conn.(s).(t).(p)),'filepath',char(p2l.conn.(s).(t)));
                icNums = [];
                for j=transpose(string(fieldnames(clustNames)))
                    if any(s==transpose(string(fieldnames(singClust.(j)))))
                        % because it is odrdered, we can add them primitvely, best practice is to look up though
                        icNums = [icNums find(j==transpose(string(fieldnames(clustNames))))];
                    end
                end
                icNums = [icNums numComp-length(muscles)+1:numComp]; % add the muscle comps becasue everybody has them :)
                numConns_boot(icNums,icNums) = numConns_boot(icNums,icNums)+1; % this concatenates the ICs in a cluster
                for j = icNums % now pull the componenets to the clusters
                    for k = icNums
                        connStruct_boot{j,k}(numConns_boot(j,k),:,:)=squeeze(EEG.CAT.Conn.dDTF08(j==icNums,k==icNums,:,:));
                        connStructS_boot{j,k}(numConns_boot(j,k),:,:)=squeeze(EEG.CAT.Conn.dDTF08(j==icNums,k==icNums,:,:));
                    end
                end
            catch
                cprintf('Yellow',s+"_"+t+"_"+p+" was not available \n")
            end
        end
        latencies = EEG.CAT.Conn.erWinCenterTimes;
%         save(p2l.studyConn+"concatConn_"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct_boot","connStructS_boot","latencies")
    end
end

%% remove baseline and create matrices for chord plot
timeRange = [-0.4 0.4];  % 0.4 sec around the event is being selected for analysis
timepts = EEG.CAT.Conn.erWinCenterTimes;
timeIdx = getindex(timepts,timeRange);
timeIdx = timeIdx(1)+1:timeIdx(2);
thetaIdx = find(EEG.CAT.Conn.freqs>=4 & EEG.CAT.Conn.freqs<=8);
alphaIdx = find(EEG.CAT.Conn.freqs>=8 & EEG.CAT.Conn.freqs<=13);
for t = trialTypes
    for p = phases
        load(p2l.studyConn+"concatConn_"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct_boot","connStructS_boot","latencies")
        load(p2l.studyConn+"concatConn_baseSub"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct","latencies")
        baseidx = find(latencies>=timewindow(1) & latencies<=timewindow(2));   % [-0.4 -0.1]
        netVals.theta=cell(numComp,numComp); netVals.alpha=cell(numComp,numComp); %netVals.all=cell(16,16);
        netVals_ave.theta=zeros(numComp,numComp); netVals_ave.alpha=zeros(numComp,numComp); %netVals_ave.all=zeros(16,16);
        for j = 1:numComp
            for k = 1:numComp
                baseVals=median(connStruct_boot{j,k}(:,:,baseidx),3);
                curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
                bootDat_theta=curr_ersp(:,thetaIdx,timeIdx);
                bootDat_alpha=curr_ersp(:,alphaIdx,timeIdx);
                for ii = 1:size(bootDat_theta,1)
                    A=bootDat_theta(ii,:,:);
                    A(squeeze(connStruct(j,k,thetaIdx,timeIdx))==0)=0; %mask using average mask
                    netVals.theta{j,k}=[netVals.theta{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,thetaIdx,timeIdx);
                    netVals_ave.theta(j,k)= mean(tmp_dat(:));
                    
                    A=bootDat_alpha(ii,:,:);
                    A(squeeze(connStruct(j,k,alphaIdx,timeIdx))==0)=0; %mask using average mask
                    netVals.alpha{j,k}=[netVals.alpha{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,alphaIdx,timeIdx);
                    netVals_ave.alpha(j,k)= mean(tmp_dat(:));
                end
            end
        end
        save(p2l.studyConn + "connMatChord" + fs +"netVals_"+groupName+"_"+t+"_"+p+"_aveTime_sbjs.mat",'netVals');
        save(p2l.studyConn + "connMatChord" + fs +"netVals_"+groupName+"_"+t+"_"+p+"_aveTime.mat",'netVals_ave');

    end
end




