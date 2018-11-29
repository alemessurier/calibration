function [ hitRate_bw, FArate_bw ] = threshold_detect_np( timeSeries,thresh_scale,delay_thresh)
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
fns=fieldnames(timeSeries.deltaF);
spikeTimes=timeSeries.spikeTimes;
deltaF=timeSeries(1).deltaF;

% initiate outputs
detect_all=[];
false_detect_all=[];

for K=1:length(fns) % loop through each movie
    fn=fns{K};
    sampRateIm=timeSeries.ImMetadata.(fn).sampRateIm;
    sampRateCA=timeSeries.ImMetadata.(fn).sampRateCA;
    spikeTrain=zeros(1,round((length(timeSeries.roiTS)/sampRateIm)*sampRateCA));
    spikeTrain(spikeTimes.(fn))=1; % make logical spike train from vector of spike times
    
    
    
    % only use movies that have cell attached recordings of equal or
    % greater length
    if ~isempty(spikeTimes.(fn)) && length(spikeTrain)/sampRateCA>=length(deltaF.(fn){1})/sampRateIm
        tmpSpTimes=spikeTimes.(fn)/sampRateCA; %convert from samples to times
        SpTimesUse=tmpSpTimes(tmpSpTimes<(length(deltaF.(fn){1})/sampRateIm-3) & tmpSpTimes>0.65); % don't use spikes that occur within 1st 0.65 sec and last 3 sec of movie
        % SpSamplesUse=spikeTimes.(fn)(tmpSpTimes<(length(deltaF.(fn){J})/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_pre=[0; diff(tmpSpTimes)];
        spacing_pre=spacing_pre(tmpSpTimes<(length(deltaF.(fn){1})/sampRateIm-3) & tmpSpTimes>0.65);
        spacing_post=[diff(tmpSpTimes); 0];
        spacing_post=spacing_post(tmpSpTimes<(length(deltaF.(fn){1})/sampRateIm-3) & tmpSpTimes>0.65);
        
        iso_inds=spacing_pre>2 & spacing_post>1; %find isolated single spikes
        empty_inds=spacing_pre>3; %find epochs with no spikes for 3 seconds
        
        iso_spikes=SpTimesUse(iso_inds); % isolated spike times in seconds
        no_spikes=SpTimesUse(empty_inds)+1; % empty spike times in seconds
        
        for J=1:length(deltaF.(fn))
            
            threshold=min(deltaF.(fn){J})+(max(deltaF.(fn){J})-min(deltaF.(fn){J}))*thresh_scale; % determine threshold by multiplying the max diff of deltaF trace by scalar, add to min.
            
            over_thresh=deltaF.(fn){J}>threshold; % frames that are over threshold
            
            detected_spikes{J}=zeros(1,length(iso_spikes)); % initiate detected spike vector
            for i=1:length(iso_spikes)
                spike=floor(iso_spikes(i)*sampRateIm); % convert spike time to movie frame #
                % if there was a threshold crossing within the allowed delay
                % periods after a spike, count the spike as detected
                if sum(over_thresh(spike:spike+delay_thresh*sampRateIm))>0
                    detected_spikes{J}(i)=1;
                end
            end
            
            false_detected{J}=zeros(1,length(no_spikes));
            
            for j=1:length(false_detected{J})
                spike=floor(no_spikes(j)*sampRateIm); % convert empty epoch time to movie frame #
                % if there was a threshold crossing within the allowed delay
                % period, count as false positive
                if sum(over_thresh(spike:spike+delay_thresh*sampRateIm))>0
                    false_detected{J}(j)=1;
                end
            end
        end
        % add results from this movie to full output vectors
        detect_all.(fn)=detected_spikes;
        false_detect_all.(fn)=false_detected;
    end
    
end

fns=fieldnames(detect_all);
for i=1:length(detect_all.(fns{1}))
    detected{i}=cellfun(@(x)detect_all.(x){i},fns,'Uni',0);
    false_detected{i}=cellfun(@(x)false_detect_all.(x){i},fns,'Uni',0);
end
false_detected_bw=cellfun(@(x)cat(2,x{:}),false_detected,'Uni',0);
detected_bw=cellfun(@(x)cat(2,x{:}),detected,'Uni',0);
hitRate_bw=cellfun(@(x)sum(x)/length(x),detected_bw,'Uni',1);
FArate_bw=cellfun(@(x)sum(x)/length(x),false_detected_bw,'Uni',1);
