function [ singleSp_template_all, burst_template,allSpikes_template_all,Fig,bootStrap_template ] = make_spikeTemplate( spikeTimes,Metadata,filtSweep,deltaF )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

singleSp_template_all=[];
burst_template_all=[];
allSpikes_template_all=[];
bootStrap_template=[];
fns=fieldnames(spikeTimes);
fns_use_inds=cellfun(@(x)round(Metadata.(x).sampRateIm)==7,fns,'Uni',1);
fns_use=fns(fns_use_inds);

for J=1:length(fns_use)
    fn=fns_use{J};
    sampRateIm=Metadata.(fn).sampRateIm;
    sampRateCA=Metadata.(fn).sampRateCA;


    
    
    
    spikeTimesAligned=spikeTimes.(fn)/sampRateCA*sampRateIm;
    
    Fig(J)=figure;
    s=subplot(2,1,1);
    plot((1:length(filtSweep.(fn)))/sampRateCA,filtSweep.(fn))
    if ~isempty(spikeTimes.(fn))
        vline(spikeTimes.(fn)/sampRateCA);
        spacing_f=[1.1; diff(spikeTimes.(fn))/sampRateCA];
        spacing_b=[0; abs(diff(flipud(spikeTimes.(fn))))/sampRateCA ];
        spacing_b=flipud(spacing_b);
        singleSp=spikeTimes.(fn)(spacing_f>1&spacing_b>1);
        burstSp=spikeTimes.(fn)(spacing_f>1);
        burstSp=burstSp(~ismember(burstSp,singleSp));
        allSpikes=spikeTimes.(fn);
        if ~isempty(singleSp)
            vline_orig(singleSp/sampRateCA,'b:')
        else
        end
        if ~isempty(burstSp)
            vline_orig(burstSp/sampRateCA,'g:')
        else
        end
    else
        singleSp=[];
        burstSp=[];
        allSpikes=[];
    end
    
    ylabel('mV');
    
    h=subplot(2,1,2);
    %         cn=cellNames{c};
    plot((1:length(deltaF.(fn).ROI1))/sampRateIm,deltaF.(fn).ROI1)
    if ~isempty(spikeTimes.(fn))
        vline(spikeTimes.(fn)/sampRateCA)
        if ~isempty(singleSp)
            vline_orig(singleSp/sampRateCA,'g:')
        else
        end
        if ~isempty(burstSp)
            vline_orig(burstSp/sampRateCA,'b:')
        end
    else
    end
    xlabel('time (sec)')
    ylabel('dF/F')
    
    
    template_int=floor(3*sampRateIm);
    singleSpUse=singleSp(singleSp<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & singleSp>0.65*sampRateCA);
    if ~isempty(singleSpUse)
        singleSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),singleSpUse,'Uni',0);
%         spike_template=cell2mat(singleSp_template);
    else
        spike_template=[];
    end
    
    burstSpUse=burstSp(burstSp<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & burstSp>0.65*sampRateCA);
    if ~isempty(burstSpUse)
        burstSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),burstSpUse,'Uni',0);
%         burst_template=cell2mat(burstSp_template);
    else
        burstSp_template=[];
    end
    
    allSpikesUse=allSpikes(allSpikes<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & allSpikes>0.65*sampRateCA);
    if ~isempty(allSpikesUse)
        allSpikes_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),allSpikesUse,'Uni',0);
%         burst_template=cell2mat(burstSp_template);
    else
        allSpikes_template=[];
    end
    
    singleSp_template_all=[singleSp_template_all; singleSp_template];
    burst_template_all=[burst_template_all; burstSp_template];
    allSpikes_template_all=[allSpikes_template_all; allSpikes_template];
    
    %% determine which spikes are detected above noise
    if ~isempty(allSpikesUse)
        [bootStrapDist,confVal,meanResponses ] = make_bootStrapDist(deltaF.(fn).ROI1,length(allSpikesUse),1:length(deltaF.(fn).ROI1),1000,sampRateIm );
        bootStrap_template=[bootStrap_template; meanResponses];
        
    end
end

linkaxes([s,h],'x')

%% find average dF/F for single spikes (spaced >1 sec apart)

sampRateIm=Metadata.(fns{1}).sampRateIm;
singleSp_template_all=horzcat(singleSp_template_all{:})';
for i=1:size(singleSp_template_all,1)
   singleSp_template_all(i,:)=singleSp_template_all(i,:)-(mean(singleSp_template_all(i,1:floor(0.5*sampRateIm))));%)/mean(singleSp_template_all(i,1:ceil(0.552*sampRateIm)));
end
figure
plot((0:size(singleSp_template_all,2)-1)/sampRateIm,singleSp_template_all)
hold on
plot((0:size(singleSp_template_all,2)-1)/sampRateIm,mean(singleSp_template_all,1),'k*-','LineWidth',2)
vline(floor(0.5*sampRateIm)/sampRateIm)
title('single spikes')
% save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','singleSp_template_all')
figure
imagesc(singleSp_template_all)
title('single spikes')
%% find average dF/F for bursts 


burst_template=horzcat(burst_template_all{:})';
for i=1:size(burst_template,1)
   burst_template(i,:)=burst_template(i,:)-(mean(burst_template(i,1:floor(0.5*sampRateIm))));%)/mean(burst_template(i,1:ceil(0.552*sampRateIm)));
end
figure
plot((0:size(burst_template,2)-1)/sampRateIm,burst_template)
hold on
plot((0:size(burst_template,2)-1)/sampRateIm,mean(burst_template,1),'k*-','LineWidth',2)
vline(floor(0.5*sampRateIm)/sampRateIm)
title('bursts')
% save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','burst_template')
figure
imagesc(burst_template)
title('bursts')
%% find average dF/F for all spikes 

allSpikes_template_all=horzcat(allSpikes_template_all{:})';
for i=1:size(allSpikes_template_all,1)
   allSpikes_template_all(i,:)=allSpikes_template_all(i,:)-(mean(allSpikes_template_all(i,1:floor(0.5*sampRateIm))));%)/mean(burst_template(i,1:ceil(0.552*sampRateIm)));
end
figure
plot((0:size(allSpikes_template_all,2)-1)/sampRateIm,allSpikes_template_all)
hold on
plot((0:size(allSpikes_template_all,2)-1)/sampRateIm,mean(allSpikes_template_all,1),'k*-','LineWidth',2)
vline(floor(0.5*sampRateIm)/sampRateIm)
title('bursts')
% save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','burst_template')
figure
imagesc(allSpikes_template_all)
title('all spikes')

% burstTimesAligned=burstSp*ephys.deltaT*Metadata.sampRate;
% burstTimesAligned=burstTimesAligned(burstTimesAligned>floor(0.5*Metadata.sampRate) & burstTimesAligned<(length(aligned_dF.ROI1)-floor(3*Metadata.sampRate)));
% singleSpikesAligned=singleSp*ephys.deltaT*Metadata.sampRate;
% singleSpikesAligned=singleSpikesAligned(singleSpikesAligned>floor(0.5*Metadata.sampRate) & singleSpikesAligned<(length(aligned_dF.ROI1)-floor(3*Metadata.sampRate)));
%
%
% [bootStrapDist,confVal,meanResponses ] = make_bootStrapDist(aligned_dF.ROI1,length(spikeTimesAligned)+1,1:length(aligned_dF.ROI1),1000,Metadata.sampRate );
% figure
% hist(bootStrapDist.MeanDF)
%
% %% calculate # of spikes (start of burst or single spike) for which dF/F >95% of mean spont. activity
%
% evokedDF_single=arrayfun(@(x)subtract_bl(aligned_dF.ROI1,x,Metadata.sampRate,0.5,3 ),singleSpikesAligned,'Uni',1);
% detected_single=evokedDF_single>confVal.meanDF;
% evokedDF_burst=arrayfun(@(x)subtract_bl(aligned_dF.ROI1,x,Metadata.sampRate,0.5,3 ),burstTimesAligned,'Uni',1);
% detected_burst=evokedDF_burst>confVal.meanDF;
% detections_single=[detections_single; detected_single];
% detections_burst=[detections_burst; detected_burst];
%
% figure
% plot(aligned_dF.ROI1)
%
% vline_orig(spikeTimesAligned)
% if ~isempty(burstTimesAligned(detected_burst))
%     vline_orig(burstTimesAligned(detected_burst),'g')
% end
%
% if ~isempty(singleSpikesAligned(detected_single))
%     vline_orig(singleSpikesAligned(detected_single),'b')
% end
%
%
% figure
% plot(mean(meanResponses,1),'k','LineWidth',1.5)
% hold on
% plot(mean(bursts_all,1),'g','LineWidth',1.5)
% plot(mean(singleSpikes_all,1),'b','LineWidth',1.5)
%
% spikeTrain_tmp=zeros(1,length(ephys.chan1.samples));
% spikeTrain_tmp(spikeTimes)=1;
% spikeTrain=[spikeTrain, spikeTrain_tmp];
%
% Traces.(fn).bursts=bursts_all;
% Traces.(fn).sampRate=Metadata.sampRate;
% Traces.(fn).single=singleSpikes_all;
% Traces.(fn).spont=meanResponses;
%

end

