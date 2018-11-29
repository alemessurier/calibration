function [ tuningData ] = make_spikeTuningPlots( caDir,whisk,stim_param,sampRateCA,spikesByWhisk )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

h1=figure;
h2=figure;
h3=figure;
set(h3,'Visible','off')
sweepNum=1;
delay=stim_param.stim.preDelay(sweepNum);
isi=stim_param.stim.isi(sweepNum);
riseFall=stim_param.stim.riseFall(sweepNum);
dur=stim_param.stim.dur(sweepNum);

stimTimes(1)=delay;
tmp=([1:stim_param.stim.numImpulses-1]*(dur+isi))+delay;
stimTimes=[stimTimes tmp]/1000;

 
for J=1:length(whisk)
    whisker=whisk{J};
    whisk_data=spikesByWhisk.(whisker);
    tmpspikes=[];
    set(0,'CurrentFigure',h1);
    h1a(J)=subplot(3,3,J);
    for i=1:size(whisk_data,1)
        [~,inds]=find(whisk_data(i,:)==1);
        tmpspikes=[tmpspikes inds];
    end
    
    edges=0:100:size(whisk_data,2);
    binCount=histc(tmpspikes,edges);
    pSp=binCount/size(whisk_data,1);
    bar(edges,pSp,'histc')
    hold on
    smoothed=smooth(pSp,5,'lowess');
     plot(edges,smoothed,'g')
      vline(stimTimes*sampRateCA,'g:');
%   smoothed=interp1(edges,smoothed,edges(1):edges(end)/size(whisk_data,2):edges(end));


    maxBin(J)=max(binCount);
    
    title(whisk{J},'FontWeight','Bold','FontSize',14)
    evoked_spikes_bs=[];
    evoked_spikes_raw=[];
    
    set(0,'CurrentFigure',h2)
    h2a(J)=subplot(3,3,J)
    plotSpikeRaster(logical(whisk_data),'PlotType','vertline');
    hold on
    vline(stimTimes*sampRateCA,'g:');
    %set(0,'CurrentFigure',h3)
    numspikes=sum(whisk_data,1);
    baseline=numspikes(1:stimTimes(1)*sampRateCA);
    
    spon_byTrial=sum(whisk_data(:,1:stimTimes(1)*sampRateCA),2);
    %     baseline=smoothed(binEdges<(stimTimes(1)-5)*sampRate);
    
    for k=1:(length(baseline)-0.07*sampRateCA-1)
        mean_sp_bin(k)=sum(baseline(k:k+0.07*sampRateCA));
    end
    blnumspikes=mean(mean_sp_bin);
    sd_blfr=std(mean_sp_bin);
    blnumspikes_1sec=blnumspikes/(0.07);
    
    stimWindowsAll=cell(1,length(stimTimes));
    for j=1:length(stimTimes)
        stimWindow=ceil(stimTimes(j)*sampRateCA):ceil((stimTimes(j)+0.07)*sampRateCA); %100 ms window
        stimWindowsAll{j}=stimWindow;
        e_spikes=sum(sum(whisk_data(:,stimWindow)));
       
        e_spikes_norm=e_spikes-blnumspikes;%+3*sd_blfr);
        evoked_spikes_bs=[evoked_spikes_bs; e_spikes_norm];
        evoked_spikes_raw=[evoked_spikes_raw; e_spikes];
        % find latency to spike for each stim
        
        binSize=edges(2)/(stim_param.settings.fskHz*1000);
        minEvoked=poissinv(.95,binSize*blnumspikes_1sec);
        sd_blfr_hz=sd_blfr/(0.07);
        
        %         minEvoked=binSize*(blfr_hz+3*sd_blfr_hz);
        bins_aboveThresh=find(smoothed(edges>stimWindow(1) & edges<stimWindow(end))>minEvoked);
        
        tmp=[ diff(bins_aboveThresh)];
        bin_inds=find(tmp==1);
        if isempty(bin_inds)
            blah=NaN;
        else
            inds_late=bins_aboveThresh(bin_inds(1));
            binEdges_thisStim=edges(edges>stimWindow(1) & edges<stimWindow(end));
            blah=binEdges_thisStim(inds_late);
        end
        

        
%                [~, blah]=find(smoothed(stimTimes(j)*sampRateCA:(stimTimes(j)+0.070)*sampRateCA)>minEvoked);
               if isempty(blah)
                   latency(j)=NaN;
               else
                   latency(j)=blah(1);
               end
    end
    
    % make vector of total number of spikes during evoked period 
    stimWindowsAll=cat(2,stimWindowsAll{:});
    e_spikes_byTrial=sum(whisk_data(:,stimWindowsAll),2);
    
    % normalize e_spikes_byTrial to same amount of time as baseline measurement
    e_spikes_byTrial=e_spikes_byTrial*(length(baseline)/length(stimWindowsAll));
    
    %         resp_latency{i}=latencies/10;
%     set(0,'CurrentFigure',h1);
%     subplot(3,3,J);
%     hold on
%       vline(latency,'g:')
    evoked_byTrial.(whisk{J})=e_spikes_byTrial;
    baseline_byTrial.(whisk{J})=spon_byTrial;
    latencies(J,:)=(latency/sampRateCA)-stimTimes;
    meanSpikeCount(J)=mean(evoked_spikes_raw);
    totalEvokedSpikes(J,:)=evoked_spikes_raw;
    pSpike(J,:)=(evoked_spikes_bs)/size(whisk_data,1);
    spontRate(J)=blnumspikes_1sec/size(whisk_data,1);
end

tuningData.latencies=latencies;
tuningData.meanSpikeCount=meanSpikeCount;
tuningData.totalEvokedSpikes=totalEvokedSpikes;
tuningData.pSpike=pSpike;
tuningData.spontRate=spontRate;
tuningData.evoked_byTrial=evoked_byTrial;
tuningData.baseline_byTrial=baseline_byTrial;

set(0,'CurrentFigure',h1)
for j=1:9;
    
    set(h1a(j),'YLim',[0,1])
    set(h1a(j),'XTick',[0:5000:stim_param.settings.samples])
    set(h1a(j),'XTickLabel',[0:5000:stim_param.settings.samples]/(sampRateCA))
    set(h1a(j),'XAxisLocation','bottom','YAxisLocation','left');
end

subplot(3,3,7)
xlabel('ms','FontWeight','Bold','FontSize',14)
ylabel('numSpikes','FontWeight','Bold','FontSize',14)
annotation('textbox','Position',[0 0.9 1 0.1],'LineStyle','none','String',stim_param.recordedFN)
linkaxes(h1a,'xy')
% saveas(h1,strcat(caDir,'psth.eps'),'epsc')

set(0,'CurrentFigure',h2)
for j=1:9;
    
    set(h2a(j),'XTick',[0:5000:stim_param.settings.samples])
    set(h2a(j),'XTickLabel',[0:5000:stim_param.settings.samples]/(sampRateCA))
    set(h2a(j),'XAxisLocation','bottom','YAxisLocation','left');
end
linkaxes(h2a,'xy')
subplot(3,3,7)
xlabel('ms','FontWeight','Bold','FontSize',14)
ylabel('trial number','FontWeight','Bold','FontSize',14)
annotation('textbox','Position',[0 0.9 1 0.1],'LineStyle','none','String',stim_param.recordedFN)
cd(caDir)
% print -depsc 'raster.eps'
% linkaxes(h1a,'xy')
saveas(h2,strcat(caDir,'raster.eps'),'epsc')

% figure
% PrWhisk=mean(pSpike,2);
% calculateCentroid(PrWhisk,0.5,1)
% colormap gray
end

