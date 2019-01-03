function [ pvals ] = permuteTest_calibration( cellPath,numReps,type )
%permutation test to determine if individual cells (from calibration expts)
%respond significantly to whisker stimulus above baseline (pre-stim)

%load results for this cell
    load([cellPath,'results.mat'],'results')
    
%find evoked frames and baseline frames indices    
fnstmp=fieldnames(results.raw.ImMetadata);
      sampRateIm=results.raw.ImMetadata.(fnstmp{1}).sampRateIm;
    int(1)=ceil(0.5*sampRateIm);
    int(2)=ceil(0.5*sampRateIm)+ceil(1*sampRateIm);

    framesEvoked=int(1):int(2);
    framesBL=1:(int(1)-1)
    
whisk=results.byWhisk.whisk;
traceByStim=results.byWhisk.traceByStim;

    
    for J=1:length(whisk)
        whiskMeans{J}=mean(traceByStim.(whisk{J})(:,framesEvoked),2);
        sponMeans{J}=mean(traceByStim.(whisk{J})(:,framesBL),2);
    end
    
    sponMeans=cat(1,sponMeans{:});
    switch type
        case 'mean'
            
            parfor J=1:length(whisk)
                sponMeansUse=datasample(sponMeans,length(whiskMeans{J}),'Replace',false);
                tmpPVal(J)=permutationTest(whiskMeans{J},sponMeansUse,numReps );
            end
            
        case 'median'
            
            parfor J=1:length(whisk)
                sponMeansUse=datasample(sponMeans,length(whiskMeans{J}),'Replace',false);
                tmpPVal(J)=permutationTest_median(whiskMeans{J},sponMeansUse,numReps );
            end
    end
    
    
    
    for J=1:length(whisk)
        pvals.(whisk{J}) = tmpPVal(J);
    end
    
end

