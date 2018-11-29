
%% set directory, files to load, etc.
caDir='J:\reduced\20150824\cell4\';
stimFiles={'20150824c4_1',...
    '20150824c4_2',...
    '20150824c4_3',...
'20150824c4_4',...
    '20150824c4_5',...
    '20150824c4_6',...
    '20150824c4_7',...
    '20150824c4_8',...
    '20150824c4_9',...
    '20150824c4_10',...
    '20150824c4_11',...
    '20150824c4_12',...
'20150824c4_13',...
'20150824c4_14',...
'20150824c4_15'...
};

imDir='J:\raw\20150824\cell4\';
imFiles={'cell4_00001_grch',...
    'cell4_00002_grch',...
    'cell4_00003_grch',...
    'cell4_00004_grch',...
    'cell4_00005_grch',...
    'cell4_00006_grch',...
    'cell4_00007_grch',...
    'cell4_00008_grch',...
    'cell4_00009_grch',...
    'cell4_00010_grch',...
    'cell4_00011_grch',...
    'cell4_00012_grch',...
    'cell4_00013_grch',...
    'cell4_00014_grch',...
    'cell4_00015_grch'...
     };


[~,stim_param]=ibtRead(strcat(caDir,stimFiles{1},'.ibt'));
whisk=fliplr({'e3' 'e2' 'e1' 'd3' 'd2' 'd1' 'c3' 'c2' 'c1'});
 ephys=FifoFileRead(strcat(caDir,stimFiles{1},'.fifo'), 1, 0);
 sampRateCA=1/ephys.deltaT;
%% analysis step one

[ spikeTimes,spikesByWhisk,deltaF,rawTimeSeries,npTimeSeries,traceByStim,Metadata,filtSweep ] = calibration_step1( caDir,stimFiles,imDir,imFiles,whisk )
%% plot tuning based on spikes

[ tuningData ] = make_spikeTuningPlots( caDir,whisk,stim_param,Metadata.cell2_00001_grch.sampRateCA,spikesByWhisk );

%% Plot tuning based on imaging data
h=figure;
sampRateIm=Metadata.(imFiles{1}).sampRateIm;
plotTuningV2single(traceByStim,h,sampRateIm)
saveas(h,strcat(caDir,'tuning.fig'),'fig')
%% Create spike templates from imaging data and spiketimes

[ singleSp_template, burst_template,allSpikes_template,figs ] = make_spikeTemplate( spikeTimes,Metadata,filtSweep,deltaF,rawTimeSeries,npTimeSeries );
mkdir(caDir,'detection_plots')
for J=1:length(figs)
    saveas(figs(J),strcat(caDir,'detection_plots\',stimFiles{J},'.fig'),'fig')
end

%% create structure of results to save 

results.spikeTuning=tuningData;
results.raw.spikeTimes=spikeTimes;
results.raw.deltaF=deltaF;
results.raw.filtSweep=filtSweep;
results.raw.ImMetadata=Metadata;
results.byWhisk.whisk=whisk;
results.byWhisk.spikesByWhisk=spikesByWhisk;
results.byWhisk.traceByStim=traceByStim;
results.template.singleSp=singleSp_template;
results.template.burst=burst_template;
reults.template.allSpikes=allSpikes_template;

save(strcat(caDir,'results.mat'),'results')

