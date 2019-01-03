function [ spikeTimes,spikesByWhisk ] = make_spikesByWhisk( caDir,stimFiles,whisk )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
spikeTrain=[];


for j=1:length(whisk)
        whisker=whisk{j};
        spikesByWhisk.(whisker)=[];
end

for K=1:length(stimFiles)
    fname=stimFiles{K};
    ephys=FifoFileRead(strcat(caDir,fname,'.fifo'), 1, 0);
    [~,stim_param]=ibtRead(strcat(caDir,fname,'.ibt'));
    
    % sort stimuli by whisker
    stimTmp=ephys.chan2.samples>1;
    stimInds=diff(stimTmp);
    stimInds=find(stimInds==1);
    stimRepInds=stimInds([0; diff(stimInds)>5000]==1);
    stimRepInds=[stimInds(1); stimRepInds];
    
    stimOrder=stim_param.stim.piezoNumber(1:length(stimRepInds));
    
    for i=1:9;
        whiskNum=i-1;
        tmpInds{i}=stimOrder==whiskNum;
    end
    
        stimTimes_whisk=cellfun(@(x)stimRepInds(x),tmpInds,'Uni',0);
    
    % find spike times
        
       bl_length=floor(stim_param.stim.preDelay(1)*stim_param.settings.fskHz);
        sampRateCA=1/ephys.deltaT;
        ephysDataRaw=ephys.chan1.samples;
        normSweeps=zeros(1,length(ephysDataRaw));
        filtSweeps=zeros(1,length(ephysDataRaw));

            filtSweeps=genButterFilter(ephysDataRaw,500,6000,4,'highpass',sampRateCA);
     
        spikes=zeros(1,length(filtSweeps));
        threshold=threshold_sweep(filtSweeps,ephys);
            overThresh=filtSweeps>threshold;
            tmp=[0; diff(overThresh)];
            spikes=tmp==1;
             [spikeTimes,~]=find(spikes==1);
 
       for i=1:9
           thisWhiskInds=stimTimes_whisk{i};
           thisWhiskInds=thisWhiskInds(thisWhiskInds<(length(spikes)-sampRateCA) & thisWhiskInds>0.5*sampRateCA);
          
           for j=1:length(thisWhiskInds)
               tmpspikes(j,:)=spikes((thisWhiskInds(j)-0.5*sampRateCA):(thisWhiskInds(j)-1+sampRateCA));
           end
               
           spikesByWhisk.(whisk{i})=[spikesByWhisk.(whisk{i}); tmpspikes];
           
           clear tmpspikes
       end
           
end

end

