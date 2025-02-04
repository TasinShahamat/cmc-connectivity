function step5p5_mats4stats2(groupName, earlyFraction)
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

%% remove the baseline and bootstrap
for t = trialTypes
    for p = phases
        load(p2l.studyConn+"concatConn_"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct_boot","connStructS_boot","latencies")
        baseidx = find(latencies>=timewindow(1) & latencies<=timewindow(2));   % [-0.4 -0.1]
        parfor j = 1:numComp
            for k = 1:numComp
                baseVals=median(connStruct_boot{j,k}(:,:,baseidx),3);
                curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
                curr_ersp=permute(curr_ersp,[2 3 1]);  % Switch dimensions. For example, 10 x 30 x 347 becomes 30 x 347 x 10

                % Bootstrap and mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',1000,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);  % Median taken across subjects
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
                curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
                connStruct(j,k,:,:)=squeeze(curr_maskedersp);

%                 curr_ersp=permute(curr_ersp,[3 1 2]); % Switch back dimensions
%                 connStruct(j,k,:,:)=squeeze(curr_ersp);
                
                
                % Calculate and subtract baseline. 
                baseVals=median(connStructS_boot{j,k}(:,:,baseidx),3);
                curr_ersp = connStructS_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
                curr_ersp=permute(curr_ersp,[2 3 1]);
                
%                 Bootstrap and mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',1000,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
                curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
                connStructS(j,k,:,:)=squeeze(curr_maskedersp);

%                 curr_ersp=permute(curr_ersp,[3 1 2]);
%                 connStructS(j,k,:,:)=squeeze(curr_ersp);
%                 
            end
        end
%         save(p2l.studyConn + fs + "unmasked" + fs +"concatConn_baseSub"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct","connStructS","latencies")
    end
end

%% save conn matrices for chord plots

alpha = 0.05;
maxVal4Plot=.00025;
timeRange = [0, 4];  % 0.4 sec around the event - full window, [-0.4, -0.1] - anticipatory, [0 0.4] - reactionary
timepts = EEG.CAT.Conn.erWinCenterTimes;
timeIndices = getindex(timepts,timeRange);
timeIndices = timeIndices(1)+1:timeIndices(2);
timepts = timepts(timeIndices);

if groupName == "young"
    EEG = pop_loadset('filename',char(f2l.conn.("PS04").(t).(p)),'filepath',char(p2l.conn.("PS04").(t)));  % not quite a random set w/ conn structure. PS04 used for young, PS08 used for old.
else
    EEG = pop_loadset('filename',char(f2l.conn.("PS08").(t).(p)),'filepath',char(p2l.conn.("PS08").(t)));
end
load(p2l.studyConn+"Conn_TimeFreqGrid_template"+groupName+"_fraction"+earlyFraction+".mat",'TimeFreqGrid_template')

for t = trialTypes
    for p = phases
        load(p2l.studyConn + fs + "unmasked" + fs +"concatConn_baseSub"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct","latencies")
        EEG.CAT.Conn.dDTF08=[];
        EEG.CAT.Conn.dDTF08=connStruct; %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
        %         baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
        EEG.CAT.curComps=1:size(EEG.CAT.Conn.dDTF08,1);
        EEG.CAT.nbchan=size(EEG.CAT.Conn.dDTF08,1);
        
        origFreqs=TimeFreqGrid_template.freqValues;
        
        
        frqbands = ["theta" "alpha" "beta"];
        
        for frqband = frqbands
            if frqband == "alpha"
                freqRange = [8 13];
            elseif frqband == "theta"
                freqRange = [4 8];
            elseif frqband == "beta"
                freqRange = [14 30];
            end
            
%             if groupName == "old"
%                 brain_indices = [1:3];
%                 muscle_indices = [4:9];
%             elseif groupName == "young"
%                 brain_indices = [1:5];
%                 muscle_indices = [6:11];
%             end
            
            netVals.(frqband)=cell(numComp,numComp); 
            netVals_ave.(frqband)=zeros(numComp,numComp);
            for i = 1:numComp
                for j = 1:numComp
                    freqs = find((origFreqs>=freqRange(1)) & (origFreqs<=freqRange(2)));
                    connection_dDTF08 = mean(squeeze(connStruct(i,j,freqs,:)), 1);                      % taking mean across the frequencies
                    connection_dDTF08 = connection_dDTF08(timeIndices);                           % Selecting only the window based on timeRange
                    
                    netVals.(frqband){j, i} = connection_dDTF08;
                    netVals_ave.(frqband)(j, i) = mean(connection_dDTF08);
                end
            end
        end
        
        save(p2l.studyConn + "connMatChord" + fs +"antc_nobase_unmaskedConnStruct_netVals_"+groupName+"_"+t+"_"+p+"_aveTime_sbjs.mat",'netVals');
        save(p2l.studyConn + "connMatChord" + fs +"antc_nobase_unmaskedConnStruct_netVals_"+groupName+"_"+t+"_"+p+"_aveTime.mat",'netVals_ave');
    end
end