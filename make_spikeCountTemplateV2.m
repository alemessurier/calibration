function [ templates,templateInds ] = make_spikeCountTemplateV2( results )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

templates=[];
templateInds.ISI_B=[];
templateInds.ISI_A=[];

fns=fieldnames(results.raw.filtSweep);
spikeTimes=results.raw.spikeTimes;
deltaF=results.raw.deltaF;

for J=1:length(fns)
    fn=fns{J};
    sampRateIm=results.raw.ImMetadata.(fn).sampRateIm;
    sampRateCA=results.raw.ImMetadata.(fn).sampRateCA;


    
    
    
    spikeTimesAligned=spikeTimes.(fn)/sampRateCA*sampRateIm;
    spikeTrain=zeros(1,length(results.raw.filtSweep.(fn)));
    spikeTrain(spikeTimes.(fn))=1;

    if ~isempty(spikeTimes.(fn))  
        tmp=diff(spikeTimes.(fn))/sampRateCA;
        tmpSpTimes=spikeTimes.(fn)(2:end-1)/sampRateCA;
        spacing_f=tmp(1:end-1);
        spacing_b=tmp(2:end);

                    template_int=floor(3*sampRateIm);

        spacing_f=spacing_f(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_b=spacing_b(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        
        templateInds.ISI_B=[templateInds.ISI_B; spacing_f];
        templateInds.ISI_A=[templateInds.ISI_A; spacing_b];

    SpTimesUse=tmpSpTimes(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
    spTemplates=arrayfun(@(x)deltaF.(fn).ROI1(floor((x-0.5)*sampRateIm):(floor((x-0.5)*sampRateIm)+template_int)),SpTimesUse,'Uni',0);
    
    templates=[templates; spTemplates];

    end
    

    
    
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

