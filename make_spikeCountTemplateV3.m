function [ templates,templateInds,bootStrap_template] = make_spikeCountTemplateV3( results )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

templates=[];
templateInds.ISI_pre=[];
templateInds.ISI_post=[];
templateInds.cumSumSp=[];
bootStrap_template=[];


fns=fieldnames(results.raw.filtSweep);
spikeTimes=results.raw.spikeTimes;
deltaF=results.raw.deltaF;

for J=1:length(fns)
    fn=fns{J};
    sampRateIm=results.raw.ImMetadata.(fn).sampRateIm;
    sampRateCA=results.raw.ImMetadata.(fn).sampRateCA;
% 
%     % find spike times by thresholding
%     filtSweeps=results.raw.filtSweep.(fn);
%     time=(1:length(filtSweeps))/sampRateCA;
%         threshold=threshold_sweep(filtSweeps,time);
%         overThresh=filtSweeps>threshold;
%         tmp=[0; diff(overThresh)];
%         spikes=tmp==1;
%         [spikeTimes,~]=find(spikes==1);
    
    
    spikeTimesAligned=spikeTimes.(fn)/sampRateCA*sampRateIm;
    spikeTrain=zeros(1,length(results.raw.filtSweep.(fn)));
    spikeTrain(spikeTimes.(fn))=1;

    if ~isempty(spikeTimes.(fn)) && length(spikeTrain)/sampRateCA>=length(deltaF.(fn).ROI1)/sampRateIm
        tmpSpTimes=spikeTimes.(fn)/sampRateCA;
        SpTimesUse=tmpSpTimes(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        SpSamplesUse=spikeTimes.(fn)(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_pre=[0; diff(tmpSpTimes)];
        spacing_pre=spacing_pre(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_post=[diff(tmpSpTimes); 0];
        spacing_post=spacing_post(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        
        if ~isempty(SpTimesUse)
            for i=1:length(SpSamplesUse)
                cumSumSp(i,:)=cumsum(spikeTrain(SpSamplesUse(i):(SpSamplesUse(i)+2*sampRateCA)));
            end
            templateInds.cumSumSp=[templateInds.cumSumSp; cumSumSp];
        end
        
        templateInds.ISI_pre=[templateInds.ISI_pre; spacing_pre];
        templateInds.ISI_post=[templateInds.ISI_post; spacing_post];
        
        template_int=floor(3*sampRateIm);
        spTemplates=arrayfun(@(x)deltaF.(fn).ROI1(floor((x-0.5)*sampRateIm):(floor((x-0.5)*sampRateIm)+template_int)),SpTimesUse,'Uni',0);
        
        templates=[templates; spTemplates];
        
        %% determine which spikes are detected above noise
        if ~isempty(SpTimesUse)
            [bootStrapDist,confVal,meanResponses ] = make_bootStrapDist(deltaF.(fn).ROI1,length(SpTimesUse),1:length(deltaF.(fn).ROI1),1000,sampRateIm );
            bootStrap_template=[bootStrap_template; meanResponses];
            
        end
        
    end
    
% results.raw.spikeTimes.(fn)=spikeTimes;
    clear cumSumSp
    
end
  
% %% find average dF/F for single spikes (spaced >1 sec apart)
% spNumbers=fieldnames(template);
% for K=1:length(spNumbers)
%    if ~isempty(template.(spNumbers{K}))
%     template.(spNumbers{K})=horzcat(template.(spNumbers{K}){:})';
%     for i=1:size(template.(spNumbers{K}),1)
%         template.(spNumbers{K})(i,:)=template.(spNumbers{K})(i,:)-(mean(template.(spNumbers{K})(i,1:floor(0.5*sampRateIm))));%)/mean(template.(spNumbers{K})(i,1:ceil(0.552*sampRateIm)));
%     end
%    end
% %     figure
% %     plot((0:size(template.(spNumbers{K}),2)-1)/sampRateIm,template.(spNumbers{K}))
%     hold on
%     plot((0:size(template.(spNumbers{K}),2)-1)/sampRateIm,mean(template.(spNumbers{K}),1),'k*-','LineWidth',2)
%     vline(floor(0.5*sampRateIm)/sampRateIm)
%     title('(spNumbers{K}) spikes')
%     % save(strcat(dir_reduced,'calibration.mat'),'VmTrace','Spikes','aligned_dF','template.(spNumbers{K})')
%     figure
    % imagesc(template.(spNumbers{K}))
    % title(strcat(spNumbers{K},'spikes'))
    %
end

