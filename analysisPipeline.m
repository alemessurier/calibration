%% set directories

dir_raw='E:\Data\raw\20150314\to_reg\';
dir_processed='E:\Data\raw\20150314\16bit_scaled\';
dir_reduced='E:\Data\reduced\20150314\';
dir_ROI_positions=['E:\Data\reduced\20150314\'];
%% register and scale raw data

ex_image=registration_1dir(dir_raw,dir_processed,'uint16');
imwrite(ex_image,strcat(dir_reduced,'ex_image.tif'),'tiff');

%% define ROIs 

[ROI_positions]=label_ROIs(dir_processed,dir_ROI_positions);
save(strcat(dir_reduced,'ROI_positions.mat'),'ROI_positions')
%% load ROIs

load(strcat(dir_reduced,'ROI_positions.mat'),'ROI_positions')
load(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries','npNormTimeSeries')
load(strcat(dir_reduced,'ROI_pos_by_cell.mat'),'ROI_pos_by_cell')

%% Extract raw fluorescence time series from ROIs

[ rawTimeSeries,npNormTimeSeries,ROI_pos_by_cell] = get_fluoTimeSeries( dir_processed, ROI_positions );
save(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries','npNormTimeSeries')
save(strcat(dir_reduced,'ROI_pos_by_cell.mat'),'ROI_pos_by_cell')
 %% Input experiment parameters
 
truncPoints = 1;           %What point to begin reading each data file (due to Pockels cell nonlinearity)
sampRate_im=7.23;             %Imaging rate
orig_recorded_frames = 723;
recordedSeconds = 100;%size(dataAll.(fns{1}).(cellNames{1}),1)/samp_rate;
badcells={};   % Cells that go off frame, etc
isi =2;                     % Inter Stimulus Interval
timeToFirstStim = 0.5;   %sec%% Clean up data
centerWhisk='d2';
fns = fieldnames(npNormTimeSeries);
% fns=fns(2:end);
[cellNames truncData] = cleanUpData(rawTimeSeries, badcells, truncPoints); % returns cellNames and truncated data without bad cells


for i=1:length(fns)
    fn=fns{i};
    for j=1:length(cellNames)
        cn=cellNames{j};
        padded=zeros(1,length(truncData.(fn).(cn))+200);
        padded(101:length(padded)-100)=truncData.(fn).(cn);
        padded(1:100)=repmat(median(truncData.(fn).(cn)),[1,100]);
        padded((end-99):end)=repmat(median(truncData.(fn).(cn)),[1,100]);
        tmpfilt=genButterFilter(padded,0.0000001,0.7,4,'butter_causal',7.23);
        filtData.(fn).(cn)=tmpfilt(101:(length(padded)-100));%-min(tmpfilt(101:(length(padded)-100)));
    end
end

 %% Employ Konnerth deltaF algorithm
    tau1=2;
    tau2=4;

[deltaF,baseline]=konnerth_deltaF_2014(filtData, tau1, tau2, sampRate_im);
% [deltaF,baseline]=calc_deltaF_AML(truncData,0.15);
truncTotal = orig_recorded_frames-length(deltaF.(fns{1}).(cellNames{1}));
%truncTotal = trunc+truncPoints;

% Determine length of frame blocks
totalFramesPerStack = size(deltaF.(fns{1}).(cellNames{1}),1);
recordedSecondsTrunc = totalFramesPerStack/sampRate_im;

%% match spike data to dF/F
cellAttached_path='E:\Data\cell_attached\20150314\';
if exist(strcat(dir_reduced,'calibration.mat'),'file')
    load(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes')
else
end

for k=1:length(fns)
    fn=fns{k};
    if ~exist(strcat(dir_reduced,'calibration.mat'),'file');
        ibtName=fn(4:end);
        ibtName(end-2)='_';
        ibt_path=strcat(cellAttached_path,ibtName,'.ibt');
        [data,param]=ibtRead(ibt_path);
        bl_length=floor(param.stim.preDelay(1)*param.settings.fskHz);
        sampRate=param.settings.fskHz;
        
        normSweeps=zeros(size(data.raw,1),size(data.raw,2));
        filtSweeps=zeros(size(data.raw,1),size(data.raw,2));
        for j=1:size(data.raw,1)
            % baseline=mean(data.raw(j,floor(size(data.raw,2)-500*sampRate):end));
            baseline=mean(data.raw(j,1:490*sampRate));
            normSweeps(j,:)=double(data.raw(j,:)-baseline);
            filtSweeps(j,:)=genButterFilter(normSweeps(j,:),500,6000,4,'highpass',10000);
        end
        
        spikes=zeros(size(filtSweeps,1),size(filtSweeps,2));
        threshold=plot_by_sweep(filtSweeps,param);
        for i=1:length(threshold)
            overThresh=filtSweeps(i,:)>threshold(i);
            tmp=[0 diff(overThresh)];
            spikes(i,:)=tmp==1;
            %     [~,inds]=find(spikes(i,:)==1);
            %     spikeTimes=[spikeTimes inds];
        end
        
        allSweeps=zeros(1,1e6);
        allSpikes=zeros(1,1e6);
        lengthSweep=size(spikes,2);
        allSweeps(1:lengthSweep)=filtSweeps(1,:);
        allSpikes(1:lengthSweep)=spikes(1,:);
        sweepTimes=param.sweepTimes-param.sweepTimes(1);
        for n=2:size(spikes,1)-1
            sweepStart=sweepTimes(n)*10000;
            allSweeps(sweepStart:sweepStart+lengthSweep-1)=filtSweeps(n,:);
            allSpikes(sweepStart:sweepStart+lengthSweep-1)=spikes(n,:);
        end
        
        VmTrace.(fn)=allSweeps;
        Spikes.(fn)=allSpikes;
        [~,spikeTimes]=find(allSpikes==1);
    else
        allSweeps=VmTrace.(fn);
        allSpikes=Spikes.(fn);
        [~,spikeTimes]=find(allSpikes==1);
    end
        
    for c=1:length(cellNames)
        cn=cellNames{c};
        tmpdF=zeros(1,orig_recorded_frames);
        tmpdF((truncTotal+1):orig_recorded_frames)=deltaF.(fn).(cn);
        aligned_dF.(fn).(cn)=tmpdF;
    end
    
    figure
    s(k)=subplot(length(cellNames)+1,1,1);
    plot((1:length(allSweeps))/10000,allSweeps)
%     vline(spikeTimes/10000)
     spacing_f=[0,diff(spikeTimes)/10000];
    spacing_b=[diff(flipud(spikeTimes))/10000,0];
    spacing_b=flipud(spacing_b);
    singleSp=spikeTimes(spacing_f>1&spacing_b>1);
    burstSp=spikeTimes(spacing_f>1);
    burstSp=burstSp(~ismember(burstSp,singleSp));
    if ~isempty(singleSp)
        vline_orig(singleSp/10000,'b:')
    else
    end
    vline_orig(burstSp/10000,'b:')
    for c=1:length(cellNames)
        h(k,c)=subplot(length(cellNames)+1,1,c+1);
        cn=cellNames{c};
        plot((1:orig_recorded_frames)/sampRate_im,aligned_dF.(fn).(cn))
%         vline(spikeTimes/10000)
        if ~isempty(singleSp)
            vline_orig(singleSp/10000)
        else
        end
            vline_orig(burstSp/10000,'b:')

        template_int=floor(2.5*sampRate_im);
        singleSp_template=arrayfun(@(x)aligned_dF.(fn).(cn)(floor((x/10000-0.552)*sampRate_im):(floor((x/10000-0.552)*sampRate_im)+template_int)),singleSp(singleSp<(1e6-(template_int*10000)/sampRate_im)),'Uni',0);
        spike_template(k).(cn)=cell2mat(singleSp_template');
        burstSp_template=arrayfun(@(x)aligned_dF.(fn).(cn)(floor((x/10000-0.552)*sampRate_im):(floor((x/10000-0.552)*sampRate_im)+template_int)),burstSp(burstSp<(1e6-(template_int*10000)/sampRate_im)),'Uni',0);
        burst_template(k).(cn)=cell2mat(burstSp_template');
    end
    linkaxes([s(k),h(k,:)],'x')
    linkaxes(h(k,:),'xy')

    % find average dF/F for single spikes (spaced >1 sec apart)
%     spacing_f=[0,diff(spikeTimes)/10000];
%     spacing_b=[0,diff(flipud(spikeTimes))/10000];
%     spacing_b=flipud(spacing_b);
%     singleSp=spikeTimes(spacing_f>1&spacing_b>1);
%     vline(singleSp)
    
    clear data
    clear param
    
end


%% find average dF/F for single spikes (spaced >1 sec apart)

rec_cell='ROI1';
tmp=arrayfun(@(x)x.(rec_cell),spike_template,'Uni',0);
singleSpikes_all=cell2mat(tmp');
for i=1:size(singleSpikes_all,1)
    singleSpikes_all(i,:)=singleSpikes_all(i,:)-(mean(singleSpikes_all(i,1:ceil(0.552*sampRate_im))));%)/mean(singleSpikes_all(i,1:ceil(0.552*sampRate_im)));
end
figure
plot((0:size(singleSpikes_all,2)-1)/7.23,singleSpikes_all)
hold on
plot((0:size(singleSpikes_all,2)-1)/7.23,mean(singleSpikes_all,1),'k*-','LineWidth',2)
vline(0.552)
title('single spikes')
% save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','singleSpikes_all')
figure
imagesc(singleSpikes_all)
title('single spikes')
%% find average dF/F for single spikes (spaced >1 sec apart)

rec_cell='ROI1';
tmp=arrayfun(@(x)x.(rec_cell),burst_template,'Uni',0);
bursts_all=cell2mat(tmp');
for i=1:size(bursts_all,1)
    bursts_all(i,:)=(bursts_all(i,:)-(mean(bursts_all(i,1:ceil(0.552*sampRate_im)))));%/mean(bursts_all(i,1:ceil(0.552*sampRate_im)));
end
figure
plot((0:size(bursts_all,2)-1)/7.23,bursts_all)
hold on
plot((0:size(bursts_all,2)-1)/7.23,mean(bursts_all,1),'k*-','LineWidth',2)
vline(0.552)
title('bursts')
% save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','singleSpikes_all')
figure
imagesc(bursts_all)
title('bursts')
save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','singleSpikes_all','bursts_all')

% [~,tmpinds]=sort(mean(bursts_all(:,1:ceil(0.552*sampRate_im)),2),'ascend')
% bursts_allSorted=bursts_all(tmpinds,:);
% figure 
% imagesc(bursts_allSorted)
% figure
% imagesc(bursts_allSorted(1:(end-4),:))
%% run analysis step 1

% lengths=cellfun(@(x)length(truncData.(x).ROI1),fns,'Uni',1);
% shortInds=lengths~=length(truncPoints:orig_recorded_frames);
% truncData=rmfield(truncData,fns{shortInds});
% fns=fns(~shortInds);

[tracesAll,stimFramesAll,cells,traceByStim,whiskers,stats]=step1(dir_reduced,...
    filtData,fns,sampRate_im,orig_recorded_frames,...
    recordedSeconds,isi,timeToFirstStim,centerWhisk,cellNames)

save(strcat(dir_reduced,'reduced_data.mat'),'tracesAll','cells','stimFramesAll','traceByStim','whiskers','stats','bootStrapResults')
