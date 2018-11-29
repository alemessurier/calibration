%%
close all
clear
clc
%% read in imaging data
imDir='J:\raw\20150825\cell2\';
imName='cell2_00001_grch';
imPath=strcat(imDir,imName,'.tif');
StimISI=3.04172;

[ImageArray, Metadata] = LoadTIFF_SI5(imPath);
[Metadata] = ReadTIFFHeader_SI5(Metadata);

%% read spike (.ibt) data

caDir='J:\reduced\20150825\';
cellName='20150825c2_1';
stimsName='20150825c2_1';

ephys=FifoFileRead(strcat(caDir,cellName,'.fifo'), 1, 0);
[spikeData,stim_param]=ibtRead(strcat(caDir,stimsName,'.ibt'));

%% Create Stimuli structure containing stimulus times, based on user-entered ISI and prestimulus delay 
PrestimDelay=stim_param.stim.preDelay(1)/1000;
nstimuli = floor(Metadata.FrameTimeInMovie(end) / StimISI) + 1;    % enough stimuli to fill the movie
Stimuli.IgorSweepNo = 0:1:nstimuli;      % Igor sweep number starts at zero.
Stimuli.Label = ones([nstimuli 1]);      % no knowledge of stimulus number here, so set all = 1.
Stimuli.Time = zeros([nstimuli 1]);
Stimuli.Time = Stimuli.IgorSweepNo * StimISI + PrestimDelay;   % time in sec from initial trigger to each stimulus

sweepNum=1;
delay=stim_param.stim.preDelay(sweepNum);
isi=stim_param.stim.isi(sweepNum);
riseFall=stim_param.stim.riseFall(sweepNum);
dur=stim_param.stim.dur(sweepNum);
stim_train=([0:stim_param.stim.numImpulses-1]*(dur+isi))/1000;



% Register Stimulus times to the movie, with no initial blackout period
[Stimuli, Metadata] = RegisterStimuliToMovie(Stimuli, Metadata, 0);
% Make a new version of the TIFF movie with strobe indicating stimulus onset times.
[ImageArrayLabeled]=LabelMovie(ImageArray,Metadata,Stimuli);

%% Define ROIs, extract fluorescence time series for each ROI

ROI_positions= define_ROIs_online(ImageArrayLabeled,[]);
[ rawTimeSeries,ROI_coord] = get_fluoTimeSeries_OL(ROI_positions,ImageArrayLabeled);

%% Get im expt parameters from tif metadata

sampRateIm=1/(Metadata.acqScanFramePeriod*Metadata.acqNumAveragedFrames);
orig_recorded_frames = Metadata.NumFrames;
recordedSeconds=orig_recorded_frames/sampRateIm;
Metadata.sampRate=sampRateIm;
%%
kurt_thresh=0;
tau1=0.2;
tau2=10;
cellNames=fieldnames(rawTimeSeries);

for J=1:length(cellNames)
    cn=cellNames{J};
    [deltaF_K.(cn) ] =  konnerth_deltaF_2014_cell(rawTimeSeries.(cn), tau1, tau2, sampRateIm);
     [ deltaF_S.(cn)] = svoboda_deltaF_simple( rawTimeSeries.(cn),kurt_thresh );

end
truncTotal = orig_recorded_frames-length(deltaF_S.(cellNames{1}));

%% match spike data to dF/F

%          ibt_path=strcat(caDir,cellName,'.ibt');
%         [data,param]=ibtRead(ibt_path);
       bl_length=floor(stim_param.stim.preDelay(1)*stim_param.settings.fskHz);
        sampRateCA=1/ephys.deltaT;
        ephysDataRaw=ephys.chan1.samples;
        normSweeps=zeros(1,length(ephysDataRaw));
        filtSweeps=zeros(1,length(ephysDataRaw));
%   
%             baseline=mean(ephysDataRaw(1:490*sampRateCA));
%             normSweeps(j,:)=double(e-baseline);
            filtSweeps=genButterFilter(ephysDataRaw,500,6000,4,'highpass',sampRateCA);
     
        spikes=zeros(1,length(filtSweeps));
        threshold=threshold_sweep(filtSweeps,ephys);
%         for i=1:length(threshold)
            overThresh=filtSweeps>threshold;
            tmp=[0; diff(overThresh)];
            spikes=tmp==1;
            %     [~,inds]=find(spikes(i,:)==1);
            %     spikeTimes=[spikeTimes inds];
     
        
      
        [spikeTimes,~]=find(spikes==1);
 
        
    for c=1:length(cellNames)
        cn=cellNames{c};
        tmpdF=zeros(1,orig_recorded_frames);
        tmpdF((truncTotal+1):orig_recorded_frames)=deltaF_S.(cn);
        aligned_dF.(cn)=tmpdF;
    end
    
    figure
     s=subplot(length(cellNames)+2,1,1);
    plot(ephys.chan1.times,filtSweeps)
     vline(spikeTimes/sampRateCA)
     spacing_f=[1.1; diff(spikeTimes)/sampRateCA];
    spacing_b=[0; abs(diff(flipud(spikeTimes)))/sampRateCA ];
    spacing_b=flipud(spacing_b);
    singleSp=spikeTimes(spacing_f>1&spacing_b>1);
    burstSp=spikeTimes(spacing_f>1);
    burstSp=burstSp(~ismember(burstSp,singleSp));
    if ~isempty(singleSp)
        vline_orig(singleSp/10000,'b:')
    else
    end
    if ~isempty(burstSp)
        vline_orig(burstSp/10000,'b:')
    else
    end
    for c=1:length(cellNames)
        h(c)=subplot(length(cellNames)+2,1,c+1);
        cn=cellNames{c};
        plot((1:orig_recorded_frames)/sampRateIm,aligned_dF.(cn))
         vline(spikeTimes/10000)
        if ~isempty(singleSp)
            vline_orig(singleSp/sampRateCA,'g:')
        else
        end
        if ~isempty(burstSp)
            vline_orig(burstSp/sampRateCA,'b:')
        end

        template_int=floor(3*sampRateIm);
        singleSpUse=singleSp(singleSp<(length(aligned_dF.(cn))-template_int)*sampRateCA/sampRateIm & singleSp>0.65*sampRateCA);
        if ~isempty(singleSpUse)
            singleSp_template=arrayfun(@(x)aligned_dF.(cn)(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),singleSpUse,'Uni',0);
            spike_template.(cn)=cell2mat(singleSp_template);
        else
            spike_template.(cn)=[];
        end
        
        burstSpUse=burstSp(burstSp<(length(aligned_dF.(cn))-template_int)*sampRateCA/sampRateIm & burstSp>0.65*sampRateCA);
        if ~isempty(burstSpUse)
            burstSp_template=arrayfun(@(x)aligned_dF.(cn)(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),burstSpUse,'Uni',0);
            burst_template.(cn)=cell2mat(burstSp_template);
        else
            burst_template.(cn)=[];
        end
    end
    h(c+1)=subplot(length(cellNames)+2,1,c+2);
    plot((1:orig_recorded_frames)/sampRateIm,rawTimeSeries.ROI1)
    vline(spikeTimes/10000)
    linkaxes([s,h],'x')
    linkaxes(h(1:c),'xy')

    % find average dF/F for single spikes (spaced >1 sec apart)
%     spacing_f=[0,diff(spikeTimes)/10000];
%     spacing_b=[0,diff(flipud(spikeTimes))/10000];
%     spacing_b=flipud(spacing_b);
%     singleSp=spikeTimes(spacing_f>1&spacing_b>1);
%     vline(singleSp)
    
  


%% find average dF/F for single spikes (spaced >1 sec apart)

rec_cell='ROI1';
tmp=arrayfun(@(x)x.(rec_cell),spike_template,'Uni',0);
singleSpikes_all=cell2mat(tmp');
for i=1:size(singleSpikes_all,1)
    singleSpikes_all(i,:)=singleSpikes_all(i,:)-(mean(singleSpikes_all(i,1:floor(0.5*sampRateIm))));%)/mean(singleSpikes_all(i,1:ceil(0.552*sampRateIm)));
end
figure
plot((0:size(singleSpikes_all,2)-1)/sampRateIm,singleSpikes_all)
hold on
plot((0:size(singleSpikes_all,2)-1)/sampRateIm,mean(singleSpikes_all,1),'k*-','LineWidth',2)
vline(floor(0.5*sampRateIm)/sampRateIm)
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
    bursts_all(i,:)=(bursts_all(i,:)-(mean(bursts_all(i,1:floor(0.5*sampRateIm)))));%/mean(bursts_all(i,1:ceil(0.552*sampRateIm)));
end
figure
plot((0:size(bursts_all,2)-1)/sampRateIm,bursts_all)
hold on
plot((0:size(bursts_all,2)-1)/sampRateIm,mean(bursts_all,1),'k*-','LineWidth',2)
vline(floor(0.5*sampRateIm)/sampRateIm)
title('bursts')
% save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','singleSpikes_all')
figure
imagesc(bursts_all)
title('bursts')
save(strcat(imDir,imName,'_calibration.mat'),'ephys','Metadata','aligned_dF','spikeTimes','singleSp','burstSp','singleSpikes_all','bursts_all')

