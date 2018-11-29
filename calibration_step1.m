function [ spikeTimes,spikesByWhisk,deltaF,rawTimeSeries,npTimeSeries,traceByStim,Metadata,filtSweep ] = calibration_step1( caDir,stimFiles,imDir,imFiles,whisk )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 if nargout>8
     error('incorrect number of outputs')
 end

if length(stimFiles)~=length(imFiles)
    error('unequal numbers of image and ephys files')
end

for a=1:length(whisk)
    whisker=whisk{a};
    spikesByWhisk.(whisker)=[];
    traceByStim.(whisker)=[];
end

for K=1:length(stimFiles)
    fname=stimFiles{K};
    fnIm=imFiles{K};
    ephys=FifoFileRead(strcat(caDir,fname,'.fifo'), 0, 0);
    [~,stim_param]=ibtRead(strcat(caDir,fname,'.ibt'));
    [ spikeTimes.(fnIm),spikesByWhisk,stimTimesWhisk,filtSweep.(fnIm) ] = spike_analysis( spikesByWhisk,whisk,ephys,stim_param );
    
    
    imPath=strcat(imDir,imFiles{K},'.tif');
    
    
    [ImageArray, Metadata.(fnIm)] = LoadTIFF_SI5(imPath);
    [Metadata.(fnIm)] = ReadTIFFHeader_SI5(Metadata.(fnIm));
    [deltaF.(fnIm),rawTimeSeries.(fnIm),npTimeSeries.(fnIm),traceByStim,truncTotal(K),Metadata.(fnIm)]=image_analysis(ImageArray,Metadata.(fnIm),stimTimesWhisk,traceByStim,whisk);
    Metadata.(fnIm).sampRateCA=1/ephys.deltaT;
   
    
end

for b=1:length(whisk)
    if ~isempty(traceByStim.(whisk{b}))
        traceByStim.(whisk{b})=horzcat(traceByStim.(whisk{b}){:})';
    end
end

    function [ spikeTimes,spikesByWhisk,stimTimesWhisk,filtSweeps ] = spike_analysis( spikesByWhisk,whisk,ephys,stim_param)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here
        
        % find spike times
        
        bl_length=floor(stim_param.stim.preDelay(1)*stim_param.settings.fskHz);
        sampRateCA=1/ephys.deltaT;
        ephysDataRaw=ephys.chan1.samples;
        normSweeps=zeros(1,length(ephysDataRaw));
        filtSweeps=zeros(1,length(ephysDataRaw));
        
        filtSweeps=genButterFilter(ephysDataRaw,500,6000,4,'highpass',sampRateCA);
        
        spikes=zeros(1,length(filtSweeps));
        time=ephys.chan2.times;
        threshold=threshold_sweep(filtSweeps,time);
        overThresh=filtSweeps>threshold;
        tmp=[0; diff(overThresh)];
        spikes=tmp==1;
        [spikeTimes,~]=find(spikes==1);
        
        % sort stimuli by whisker
        stimTmp=ephys.chan2.samples>1;
        stimInds=diff(stimTmp);
        stimInds=find(stimInds==1);
        stimRepInds=stimInds([0; diff(stimInds)>5000]==1);
        stimRepInds=[stimInds(1); stimRepInds];
        
        if length(stim_param.stim.piezoNumber)<length(stimRepInds) || length(stimRepInds)<length(whisk)
            for i=1:9
                stimTimesWhisk.(whisk{i})=[];
            end
            sprintf('%s',['stim file ',stim_param.fn,' does not have the correct number of stimulus labels.'])
        else
            
            
            stimOrder=stim_param.stim.piezoNumber(1:length(stimRepInds));
            
            for i=1:9;
                whiskNum=i-1;
                tmpInds{i}=stimOrder==whiskNum;
            end
            
            stimFrames_whisk=cellfun(@(x)stimRepInds(x),tmpInds,'Uni',0);
            
            
            
            % make spikesByWhisk
            for i=1:9
                thisWhiskInds=stimFrames_whisk{i};
                thisWhiskInds=thisWhiskInds(thisWhiskInds<(length(spikes)-sampRateCA) & thisWhiskInds>0.5*sampRateCA);
                
                if ~isempty(thisWhiskInds)
                    for j=1:length(thisWhiskInds)
                        tmpspikes(j,:)=spikes((thisWhiskInds(j)-0.5*sampRateCA):(thisWhiskInds(j)-1+sampRateCA));
                    end
                    
                    
                    spikesByWhisk.(whisk{i})=[spikesByWhisk.(whisk{i}); tmpspikes];
                    stimTimesWhisk.(whisk{i})=stimFrames_whisk{i}/sampRateCA; %stim times in seconds for aligning to imaging data
                else
                    stimTimesWhisk.(whisk{i})=[];
                end
            end
            
        end
        
        
        
        
        
    end


    function [deltaF,rawTimeSeries,npTimeSeries,traceByStim,truncTotal,Metadata]=image_analysis(ImageArray,Metadata,stimTimesWhisk,traceByStim,whisk)
        % Register Stimulus times to the movie, with no initial blackout period
        %         [Stimuli, Metadata] = RegisterStimuliToMovie(Stimuli, Metadata, 0);
        %         % Make a new version of the TIFF movie with strobe indicating stimulus onset times.
        %         [ImageArrayLabeled]=LabelMovie(ImageArray,Metadata,Stimuli);
        %
        %% Define ROIs, extract fluorescence time series for each ROI
        
        ROI_positions= define_ROIs_online(ImageArray,[]);
        [ rawTimeSeries,npTimeSeries,ROI_coord] = get_fluoTimeSeries_OL(ROI_positions,ImageArray);
        
        %% Get im expt parameters from tif metadata
        
        sampRateIm=1/(Metadata.acqScanFramePeriod*Metadata.acqNumAveragedFrames);
        orig_recorded_frames = Metadata.NumFrames;
        recordedSeconds=orig_recorded_frames/sampRateIm;
        Metadata.sampRateIm=sampRateIm;
        %%
        kurt_thresh=0;
        tau1=0.2;
        tau2=10;
        cellNames=fieldnames(rawTimeSeries);
        
        for J=1:length(cellNames)
            cn=cellNames{J};
            
            [ deltaF.(cn)] = svoboda_deltaF_simple( rawTimeSeries.(cn),kurt_thresh );
            
        end
        truncTotal = orig_recorded_frames-length(deltaF.(cellNames{1}));
        
        if round(sampRateIm)==7
            for i=1:length(whisk);
                if ~isempty(stimTimesWhisk.(whisk{i}));
                    stimTimesThisWhisk=stimTimesWhisk.(whisk{i});
                    stimFramesWhiskIm=stimTimesThisWhisk*Metadata.sampRateIm;
                    SFWI_use=stimFramesWhiskIm(stimFramesWhiskIm>0.5*Metadata.sampRateIm & stimFramesWhiskIm<(length(deltaF.ROI1)-3*Metadata.sampRateIm));
                    tmpTraceByStim=arrayfun(@(x)deltaF.ROI1...
                        ((x-floor(0.5*Metadata.sampRateIm)):(x+floor(2.5*Metadata.sampRateIm)))...
                        -mean(deltaF.ROI1((x-floor(0.5*Metadata.sampRateIm)):x)),SFWI_use,'Uni',0);
                    traceByStim.(whisk{i})=[traceByStim.(whisk{i}); tmpTraceByStim];
                end
            end
        else
        end
    end


   




end



