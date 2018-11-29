function [ template ] = make_spikeCountTemplate( results,timeInt )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 template.single=[];
       template.double=[];
       template.triple=[];
       template.four=[];
       template.five=[];


fns=fieldnames(results.raw.filtSweep);
spikeTimes=results.raw.spikeTimes;
deltaF=results.raw.deltaF;

for J=1:length(fns)
    fn=fns{J};
    sampRateIm=results.raw.ImMetadata.(fn).sampRateIm;
    sampRateCA=results.raw.ImMetadata.(fn).sampRateCA;

clear numSpikes
    
    
    
    spikeTimesAligned=spikeTimes.(fn)/sampRateCA*sampRateIm;
    spikeTrain=zeros(1,length(results.raw.filtSweep.(fn)));
    spikeTrain(spikeTimes.(fn))=1;
%     
%     Fig(J)=figure;
%     s=subplot(2,1,1);
%     plot((1:length(filtSweep.(fn)))/sampRateCA,filtSweep.(fn))
%     
    if ~isempty(spikeTimes.(fn))  
       spacing_f=[1.1; diff(spikeTimes.(fn))/sampRateCA];
        firstSp=spikeTimes.(fn)(spacing_f>1);
        firstSp=firstSp(firstSp<(length(spikeTrain)-1*sampRateCA));
        if ~isempty(firstSp)
            for i=1:length(firstSp)
                numSpikes(i)=sum(spikeTrain(firstSp(i):(firstSp(i)+timeInt*sampRateCA)));
            end
            
            
            singleSpTimes=firstSp(numSpikes==1);
            doubleSpTimes=firstSp(numSpikes==2);
            tripleSpTimes=firstSp(numSpikes==3);
            fourSpTimes=firstSp(numSpikes==4);
            fiveSpTimes=firstSp(numSpikes>=5);
            
            
            allSpikes=spikeTimes.(fn);
        end
    else
         singleSpTimes=[];
       doubleSpTimes=[];
       tripleSpTimes=[];
       fourSpTimes=[];
       fiveSpTimes=[];
    end
    
%     ylabel('mV');
%     
%     h=subplot(2,1,2);
%     %         cn=cellNames{c};
%     plot((1:length(deltaF.(fn).ROI1))/sampRateIm,deltaF.(fn).ROI1)
%     if ~isempty(spikeTimes.(fn))
%         vline(spikeTimes.(fn)/sampRateCA)
%         if ~isempty(singleSp)
%             vline_orig(singleSp/sampRateCA,'g:')
%         else
%         end
%         if ~isempty(burstSp)
%             vline_orig(burstSp/sampRateCA,'b:')
%         end
%     else
%     end
%     xlabel('time (sec)')
%     ylabel('dF/F')
    
    
    template_int=floor(3*sampRateIm);
    singleSpUse=singleSpTimes(singleSpTimes<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & singleSpTimes>0.65*sampRateCA);
    if ~isempty(singleSpUse)
        singleSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),singleSpUse,'Uni',0);
%         spike_template=cell2mat(singleSp_template);
    else
        singleSp_template=[];
    end
    
    doubleSpUse=doubleSpTimes(doubleSpTimes<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & doubleSpTimes>0.65*sampRateCA);
    if ~isempty(doubleSpUse)
        doubleSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),doubleSpUse,'Uni',0);
%         burst_template=cell2mat(burstSp_template);
    else
        doubleSp_template=[];
    end
    
     tripleSpUse=tripleSpTimes(tripleSpTimes<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & tripleSpTimes>0.65*sampRateCA);
    if ~isempty(tripleSpUse)
        tripleSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),tripleSpUse,'Uni',0);
%         burst_template=cell2mat(burstSp_template);
    else
        tripleSp_template=[];
    end
    
     fourSpUse=fourSpTimes(fourSpTimes<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & fourSpTimes>0.65*sampRateCA);
    if ~isempty(fourSpUse)
        fourSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),fourSpUse,'Uni',0);
%         burst_template=cell2mat(burstSp_template);
    else
        fourSp_template=[];
    end
    
    fiveSpUse=fiveSpTimes(fiveSpTimes<(length(deltaF.(fn).ROI1)-template_int)*sampRateCA/sampRateIm & fiveSpTimes>0.65*sampRateCA);
    if ~isempty(fiveSpUse)
        fiveSp_template=arrayfun(@(x)deltaF.(fn).ROI1(floor((x/sampRateCA-0.5)*sampRateIm):(floor((x/sampRateCA-0.5)*sampRateIm)+template_int)),fiveSpUse,'Uni',0);
%         burst_template=cell2mat(burstSp_template);
    else
        fiveSp_template=[];
    end
    
 template.single=[template.single; singleSp_template];
       template.double=[template.double; doubleSp_template];
       template.triple=[template.triple; tripleSp_template];
       template.four=[template.four; fourSp_template];
       template.five=[template.five; fiveSp_template];
    
    
end
  
%% find average dF/F for single spikes (spaced >1 sec apart)
spNumbers=fieldnames(template);
for K=1:length(spNumbers)
   if ~isempty(template.(spNumbers{K}))
    template.(spNumbers{K})=horzcat(template.(spNumbers{K}){:})';
    for i=1:size(template.(spNumbers{K}),1)
        template.(spNumbers{K})(i,:)=template.(spNumbers{K})(i,:)-(mean(template.(spNumbers{K})(i,1:floor(0.5*sampRateIm))));%)/mean(template.(spNumbers{K})(i,1:ceil(0.552*sampRateIm)));
    end
   end
%     figure
%     plot((0:size(template.(spNumbers{K}),2)-1)/sampRateIm,template.(spNumbers{K}))
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

