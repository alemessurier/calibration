function [ detect_all, false_detect_all,misses_all,correct_rej_all ] = threshold_detect_event( results,thresh_scale,delay_thresh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fns=fieldnames(results.raw.spikeTimes);
spikeTimes=results.raw.spikeTimes;
deltaF=results.raw.deltaF;

detect_all=[];
false_detect_all=[];
misses_all=[];
correct_rej_all=[];

for K=1:length(fns)
    fn=fns{K};
    sampRateIm=results.raw.ImMetadata.(fn).sampRateIm;
    sampRateCA=results.raw.ImMetadata.(fn).sampRateCA;
   
    
    threshold=median(deltaF.(fn).ROI1)+(max(deltaF.(fn).ROI1)-min(deltaF.(fn).ROI1))*thresh_scale;
    over_thresh=deltaF.(fn).ROI1>threshold;
    thresh_cross=find(diff(over_thresh)==1);
    
   
    spTimesAct=spikeTimes.(fn)/sampRateCA;
    
    
    detected_spikes=zeros(1,length(thresh_cross));
    false_detected=zeros(1,length(thresh_cross));
    for i=1:length(thresh_cross)
        tc_time=thresh_cross(i)/sampRateIm;
        spike=sum(spTimesAct>(tc_time-delay_thresh) & spTimesAct<tc_time);
          detected_spikes(i)=spike;
        if spike==0
          false_detected(i)=1;
        end
        
    end
    
     delayFrames=floor(delay_thresh*sampRateIm);
     spTimesAct=spTimesAct(floor(spTimesAct*sampRateIm)>0 & spTimesAct*sampRateIm<(length(over_thresh)-delayFrames));
    misses=zeros(1,length(spTimesAct));
   
    for j=1:length(spTimesAct)
        spike=floor(spTimesAct(j)*sampRateIm);
        if sum(over_thresh(spike:(spike+delayFrames)))==0
            misses(j)=1;
        end
    end
    
    correct_reject=zeros(1,length(thresh_cross));
    below_thresh=~over_thresh;
    inds_below_thresh=(1:length(deltaF.(fn).ROI1));
    inds_below_thresh=inds_below_thresh(below_thresh);
    bt_samp=inds_below_thresh(randi(length(inds_below_thresh),1,length(thresh_cross)));
    
    for j=1:length(bt_samp)
        bt_time=floor(bt_samp(j)/sampRateIm);
        spike=sum(spTimesAct>(bt_time-delay_thresh) & spTimesAct<bt_time);
        if spike==0
            correct_reject(j)=1;
        end
    end
    
    detect_all=[detect_all detected_spikes];
    false_detect_all=[false_detect_all false_detected];
    misses_all=[misses_all misses];
    correct_rej_all=[correct_rej_all correct_reject];
end

end

