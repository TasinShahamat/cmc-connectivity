function step5_group_connectvity_plot(subjs, groupName, earlyFraction)
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
if ~exist('groupName','var') || isempty(groupName), groupName = "old"; else, groupName = string(groupName); end
if ~exist('earlyFraction','var') || isempty(earlyFraction), earlyFraction = 1; end % fraction of pertubation to be flagged as early (or late)
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
trialTypes = ["LEI", "LME", "REI", "RME"];
if earlyFraction < .5
    phases = ["early_epoch_CPT" "early_epoch_PPT" "late_epoch_CPT" ...
        "late_epoch_PPT"];
else
    phases = ["early_epoch_CPT" "early_epoch_PPT"]; % just running connectivty for perturbed and catch 
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
p2l.eeglab = p2l.git + fs + "eeglab2019_0" + fs;
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
        save(p2l.studyConn+"concatConn_"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct_boot","connStructS_boot","latencies")
    end
end

%% remove the baseline and bootstrap
for t = trialTypes
    for p = phases
        load(p2l.studyConn+"concatConn_"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct_boot","connStructS_boot","latencies")
        baseidx = find(latencies>=timewindow(1) & latencies<=timewindow(2));
        for j = 1:numComp
            for k = 1:numComp
                % Calculate and subtract baseline
                baseVals=median(connStruct_boot{j,k}(:,:,baseidx),3);
                curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
                curr_ersp=permute(curr_ersp,[2 3 1]);
                
                % Bootstrap and mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',200,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
                curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
                connStruct(j,k,:,:)=squeeze(curr_maskedersp);
                
                % Calculate and subtract baseline
                baseVals=median(connStructS_boot{j,k}(:,:,baseidx),3);
                curr_ersp = connStructS_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
                curr_ersp=permute(curr_ersp,[2 3 1]);
                
                % Bootstrap and mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',200,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
                curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
                connStructS(j,k,:,:)=squeeze(curr_maskedersp);
                
            end
        end
        save(p2l.studyConn+"concatConn_baseSub"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct","connStructS","latencies")
    end
end

%% plot
% We need to benfent from the EEG.CAT strucutre to plot the connectivity
% matrix. ;)
maxVal4Plot=.00025;
EEG = pop_loadset('filename',char(f2l.conn.("PS08").(t).(p)),'filepath',char(p2l.conn.("PS08").(t)));  % not quite a random set w/ conn structure
load(p2l.studyConn+"Conn_TimeFreqGrid_template"+groupName+"_fraction"+earlyFraction+".mat",'TimeFreqGrid_template')
for t = trialTypes
    for p = phases
        load(p2l.studyConn+"concatConn_baseSub"+groupName+"_"+t+"_"+p+"_fraction"+earlyFraction+".mat","connStruct","latencies")
        EEG.CAT.Conn.dDTF08=[];
        EEG.CAT.Conn.dDTF08=connStruct; %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
%         baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
        EEG.CAT.curComps=1:size(EEG.CAT.Conn.dDTF08,1);
        EEG.CAT.nbchan=size(EEG.CAT.Conn.dDTF08,1);
%         EEGtemp=EEG;
%         EEGtemp = pop_vis_TimeFreqGrid(EEGtemp,'nogui','topoplot','');
%         close(gcf)
        EEG.CAT.configs.vis_TimeFreqGrid=TimeFreqGrid_template;  % another cool hacking :)
        EEG.CAT.configs.vis_TimeFreqGrid.nodelabels=[cellstr(bNames) cellstr(muscles)];
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.clim=[-maxVal4Plot maxVal4Plot];
% %
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout=[];
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.arg_direct=0;
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.triu='dDTF08';
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.ut_clim=[-maxVal4Plot maxVal4Plot];
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.tril='dDTF08';
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.lt_clim=[-maxVal4Plot maxVal4Plot];
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.diag='none';
        maxValMeasS=3; %1.5;
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.d_clim=[-maxValMeasS maxValMeasS];
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.clim=[];
        EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.arg_selection='Partial';
        EEG.CAT.configs.vis_TimeFreqGrid.timeRange=[];
        EEG.CAT.configs.vis_TimeFreqGrid.smooth = 0; % changing to one may change the figure altogether 
        EEG.CAT.configs.vis_TimeFreqGrid.colormap = 'parula(300)';
        EEG.CAT.configs.vis_TimeFreqGrid.freqscale='linear';%'linear';%'log'; % %Note: log scale doesn't seem to alter plotting result
        EEG.CAT.configs.vis_TimeFreqGrid.foilines=[];%log([4 8 13 30]); %
        EEG.CAT.configs.vis_TimeFreqGrid.foilinecolor=[0 0 0]; %[1 0 0; 1 0 0; 0 0 1; 0 0 1];
        origFreqs=TimeFreqGrid_template.freqValues;
        EEG.CAT.configs.vis_TimeFreqGrid.freqValues=origFreqs(origFreqs>=4);
        EEG.CAT.configs.vis_TimeFreqGrid.events{1}={0 'k' ':' 2};
        if strcmp(cond(1),'a')
            EEG.CAT.configs.vis_TimeFreqGrid.events{2}={0.5 'k' ':' 2};
        else
            EEG.CAT.configs.vis_TimeFreqGrid.events{2}={1 'k' ':' 2};
        end
        EEG.CAT.configs.vis_TimeFreqGrid.timeRange=[-0.4 0.4];
        EEG = pop_vis_TimeFreqGrid(EEG,'nogui','topoplot','');
        set(gcf,"Name",t+"_"+p+"_masked_"+groupName+"_fraction"+earlyFraction,'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
        saveas(gcf, p2l.studyConn+"figs/" + t+"_"+p+"_masked_"+groupName+"_fraction"+earlyFraction+".fig")
        print(gcf, p2l.studyConn+"figs/" + t+"_"+p+"_masked_"+groupName+"_fraction"+earlyFraction+ ".pdf", "-dpdf", "-r300", "-painters", "-bestfit");
        print(gcf, p2l.studyConn+"figs/" + t+"_"+p+"_masked_"+groupName+"_fraction"+earlyFraction+".png", "-dpng", "-r300", "-painters");
    end
end
































