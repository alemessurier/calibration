function [ detect_all, false_detect_all ] = threshold_detect( results,thresh_scale,delay_thresh)
% THRESHOLD_DETECT thresholds deltaF/F data and determines correct spike
% detections and false positives by comparison w/simultaneous spike train
% data. 
%
% INPUTS
%           'results'           structure of reduced data from one cell
%                               containing cell-attached and imaging data
%
%           'thresh_scale'      scalar value used to calculate threshold
%
%           'delay_thresh'      time period following spike in which threshold
%                               crossing will be counted as detecting that spike
%
% OUTPUTS
%           'detect_all'        binary vector that is the same length as number
%                               of spikes. 1=detected spike, 0=miss
%               
%           'false_detect_all'  binary vector that is the same length as number
%                               of epochs with no spikes. 1=threhold crossing
%                                (false positive), 0=correct rejection

% read in variables from results structure
fns=fieldnames(results.raw.spikeTimes); 
spikeTimes=results.raw.spikeTimes;
deltaF=results.raw.deltaF;

% initiate outputs
detect_all=[]; 
false_detect_all=[];

for K=1:length(fns) % loop through each movie
    fn=fns{K};
    sampRateIm=results.raw.ImMetadata.(fn).sampRateIm;
    sampRateCA=results.raw.ImMetadata.(fn).sampRateCA;
    spikeTrain=zeros(1,length(results.raw.filtSweep.(fn)));
    spikeTrain(spikeTimes.(fn))=1; % make logical spike train from vector of spike times
    
    threshold=min(deltaF.(fn).ROI1)+(max(deltaF.(fn).ROI1)-min(deltaF.(fn).ROI1))*thresh_scale; % determine threshold by multiplying the max diff of deltaF trace by scalar, add to min.
    
    % only use movies that have cell attached recordings of equal or
    % greater length
    if ~isempty(spikeTimes.(fn)) && length(spikeTrain)/sampRateCA>=length(deltaF.(fn).ROI1)/sampRateIm 
        tmpSpTimes=spikeTimes.(fn)/sampRateCA; %convert from samples to times
        SpTimesUse=tmpSpTimes(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65); % don't use spikes that occur within 1st 0.65 sec and last 3 sec of movie
       % SpSamplesUse=spikeTimes.(fn)(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_pre=[0; diff(tmpSpTimes)];
        spacing_pre=spacing_pre(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_post=[diff(tmpSpTimes); 0];
        spacing_post=spacing_post(tmpSpTimes<(length(deltaF.(fn).ROI1)/sampRateIm-3) & tmpSpTimes>0.65);
        
        iso_inds=spacing_pre>2 & spacing_post>1; %find isolated single spikes
        empty_inds=spacing_pre>3; %find epochs with no spikes for 3 seconds
        
        iso_spikes=SpTimesUse(iso_inds); % isolated spike times in seconds
        no_spikes=SpTimesUse(empty_inds)+2; % empty spike times in seconds
        
        over_thresh=deltaF.(fn).ROI1>threshold; % frames that are over threshold
        
        detected_spikes=zeros(1,length(iso_spikes)); % initiate detected spike vector
        for i=1:length(iso_spikes)
            spike=floor(iso_spikes(i)*sampRateIm); % convert spike time to movie frame #
            % if there was a threshold crossing within the allowed delay
            % periods after a spike, count the spike as detected
            if sum(over_thresh(spike:spike+delay_thresh*sampRateIm))>0
                detected_spikes(i)=1;
            end
        end
        
        false_detected=zeros(1,length(no_spikes));
        
        for j=1:length(false_detected)
             spike=floor(no_spikes(j)*sampRateIm); % convert empty epoch time to movie frame #
            % if there was a threshold crossing within the allowed delay
            % period, count as false positive
            if sum(over_thresh(spike:spike+delay_thresh*sampRateIm))>0
                false_detected(j)=1;
            end
        end
        
        % add results from this movie to full output vectors
        detect_all=[detect_all detected_spikes];
        false_detect_all=[false_detect_all false_detected];
    end
            
end

