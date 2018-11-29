%% compare tuning curve widths measured w/spikes vs. deltaF/F

 cellPaths={'J:\reduced\20150817\cell2\',...
    'J:\reduced\20150819\cell2\',...
    'J:\reduced\20150819\cell3\',...
    'J:\reduced\20150824\cell1\',...
    'J:\reduced\20150824\cell2\',...
    'J:\reduced\20150824\cell4\',...
    'J:\reduced\20150825\cell2\',...
    'J:\reduced\20150825\cell1\',...
    'J:\reduced\20150827\cell1\',...
    'J:\reduced\20150827\cell2\',...
    'J:\reduced\20150827\cell3\'...
    };

% int(1)=ceil(0.5*sampRateIm);
% int(2)=ceil(0.5*sampRateIm)+round(0.5*sampRateIm);
colormap gray
% areaSp=zeros(1,9);
% areaDF=zeros(1,9);
RF_proj=zeros(length(cellPaths),1);
for K=1:length(cellPaths)
    %load in reduced data for this cell
    load(strcat(cellPaths{K},'results.mat'));
    cd(cellPaths{K});
    
    %load in stimulus files
    name=dir('analysisPipeline_*');
    filecont=fileread(strcat(cellPaths{K},name.name));
expr = '[^\n]*stimFiles=[^\n]*,...*};';
dp_string = regexp(filecont,expr,'match');
eval(dp_string{:});
    
% get stimulus parameters
[~,stim_param]=ibtRead(strcat(cellPaths{K},stimFiles{1},'.ibt'));

    figure;
%     subplot(1,2,1)
    RF_heatmap(results.byWhisk.traceByStim)
    title(cellPaths{K})    
fnstmp=fieldnames(results.raw.ImMetadata);
       sampRateCA=results.raw.ImMetadata.(fnstmp{1}).sampRateCA;
    sampRateIm=results.raw.ImMetadata.(fnstmp{1}).sampRateIm;

    whisk=results.byWhisk.whisk;   
    spikeTuning=make_spikeTuningPlots( cellPaths{K},whisk,stim_param,sampRateCA,results.byWhisk.spikesByWhisk );
  
    int(1)=ceil(0.5*sampRateIm);
    int(2)=ceil(0.5*sampRateIm)+ceil(1*sampRateIm);

   
    responseVecSp=mean(results.spikeTuning.pSpike,2);
    
    responseVec_dF=cellfun(@(x)mean(mean(results.byWhisk.traceByStim.(x)(:,int(1):int(2)))),whisk,'Uni',1);
    
    % plot receptive field as a heatmap
   
    responseMapSp=reshape(responseVecSp,[3 3])';
    responseMapDF=reshape(responseVec_dF,[3 3])';
figure;
subplot(1,2,1)
    imagesc(responseMapSp)
colormap gray
set(gca,'XTick',[1:3]);
        set(gca,'XTickLabel',{'1' '2' '3'},'FontWeight','bold');
        set(gca,'YTick',[1:3]);
        set(gca,'YTickLabel',{'C' 'D' 'E'},'FontWeight','bold');
        set(gca,'LineWidth',2);
axis square
title('spikes')

subplot(1,2,2)
    imagesc(responseMapDF)
colormap gray
set(gca,'XTick',[1:3]);
        set(gca,'XTickLabel',{'1' '2' '3'},'FontWeight','bold');
        set(gca,'YTick',[1:3]);
        set(gca,'YTickLabel',{'C' 'D' 'E'},'FontWeight','bold');
        set(gca,'LineWidth',2);
axis square
title('dF/F')

    
    RF_proj(K)=acos(dot(responseVecSp,responseVec_dF)/(norm(responseVecSp)*norm(responseVec_dF)))*(180/pi);

    % update results file
    save([cellPaths{K},'results_old.mat'],'results')
    results.spikeTuning=spikeTuning;
    save([cellPaths{K},'results.mat'],'results')
end

%% make plots comparing BW/tuning from spike measurements and average dF/F

RFs_spike=zeros(length(cellPaths),9);
RFs_dF=RFs_spike;

for K=1:length(cellPaths)
    %load in reduced data for this cell
    load(strcat(cellPaths{K},'results.mat'));
    
    % calculate evoked frames for dF/F
    fnstmp=fieldnames(results.raw.ImMetadata);
    sampRateIm=results.raw.ImMetadata.(fnstmp{1}).sampRateIm;
    int(1)=ceil(0.5*sampRateIm);
    int(2)=ceil(0.5*sampRateIm)+ceil(1*sampRateIm);

    whisk=results.byWhisk.whisk;   
    
    responseVecSp=mean(results.spikeTuning.pSpike,2);
    semSP=std(results.spikeTuning.pSpike,[],2);
    responseVec_dF=cellfun(@(x)mean(mean(results.byWhisk.traceByStim.(x)(:,int(1):int(2)))),whisk,'Uni',1);
    
    RFs_spike(K,:)=responseVecSp';
    RFs_dF(K,:)=responseVec_dF;
end

%Zscore each RF - spikes
std_sp=std(RFs_spike,[],2);
std_sp=repmat(std_sp,1,9);
means_sp=mean(RFs_spike,2);
means_sp=repmat(means_sp,1,9);
tmp_sp=RFs_spike-means_sp;
RFs_spike_Z=tmp_sp./std_sp

%Zscore each RF - df
std_dF=std(RFs_dF,[],2);
std_dF=repmat(std_dF,1,9);
means_dF=mean(RFs_dF,2);
means_dF=repmat(means_dF,1,9);
tmp_dF=RFs_dF-means_dF;
RFs_dF_Z=tmp_dF./std_dF;

% find BWs
[~,BWs_sp]=max(RFs_spike_Z,[],2);
[~,BWs_dF]=max(RFs_dF_Z,[],2);
figure; hold on
plot(BWs_sp,BWs_dF,'ko')
tmp=gca;
tmp.XLim=[1 9];
tmp.YLim=[1 9];

BWs_mat=zeros(9,9,length(BWs_dF));
for i=1:length(BWs_dF)
    BWs_mat(BWs_sp(i),BWs_dF(i),i)=1;
end
BWs_mat=sum(BWs_mat,3);
figure; imagesc(BWs_mat);


% find angle between spike and dF RFs
for K=1:size(RFs_spike_Z,1)
    spRF=RFs_spike_Z(K,:);
    dFRF=RFs_dF_Z(K,:);
    RF_proj(K)=acos(dot(spRF,dFRF)/(norm(spRF)*norm(dFRF)))*(180/pi);
    
    % make shuffled distribution of angles
    for t=1:10000
        order=randperm(9);
        spRF=RFs_spike_Z(K,order);
        order=randperm(9);
        dFRF=RFs_dF_Z(K,order);
    
        proj_shuff(t)=acos(dot(spRF,dFRF)/(norm(spRF)*norm(dFRF)))*(180/pi);
    end
    
    tmpidx=proj_shuff>180;
    
    mean_shuff_proj(K)=mean(proj_shuff);
    std_shuff_proj(K)=std(proj_shuff);
    
    % find 5th and 95th percentile of shuffled distribution
    proj_shuff=sort(proj_shuff,'ascend');
    idx_5th=ceil(0.05*length(proj_shuff));
    shuff_5thP(K)=proj_shuff(idx_5th);
    idx_95th=ceil(0.95*length(proj_shuff));
    shuff_95thP(K)=proj_shuff(idx_95th);
    
end
shuff_95thP=shuff_95thP-mean_shuff_proj;
figure; hold on
plot(1:length(BWs_dF),RF_proj,'b*')
errorbar(1:length(BWs_dF),mean_shuff_proj,shuff_5thP,shuff_95thP,'k.')

%% 2d tuning plots
figure; hold on
for K=1:length(RFs_dF_Z)
    h(K)=subplot(6,2,K)
    hold on
    plot(1:9,RFs_dF_Z(K,:),'k.-')
    plot(1:9,RFs_spike_Z(K,:),'b.-')
    ylabel(num2str(K));
end
tmp=gca;
    tmp.XTick=1:9;
    tmp.XTickLabel=whisk
    ylabel('Z-scored response')
    legend('dF/F','spikes')
    linkaxes(h,'x')

%% 2d tuning plots (spikes, not Z-scored)
figure; hold on
for K=1:length(RFs_dF_Z)
    h(K)=subplot(6,2,K)
    hold on
    plot(1:9,RFs_spike(K,:),'b.-')
    ylabel(num2str(K));
end
tmp=gca;
    tmp.XTick=1:9;
    tmp.XTickLabel=whisk
    ylabel('Z-scored response')
    legend('dF/F','spikes')
    linkaxes(h,'x')

%% tuning idx
min_sp=min(RFs_spike_Z,[],2);
min_sp=repmat(min_sp,1,9);
RFs_sp=RFs_spike_Z+abs(min_sp);

min_dF=min(RFs_dF_Z,[],2);
min_dF=repmat(min_dF,1,9);
RFs_dF=RFs_dF_Z+abs(min_dF);

for i=1:size(RFs_dF,1)
    TI_sp(i)=max(RFs_sp(i,:))/norm(RFs_sp(i,:));
    TI_dF(i)=max(RFs_dF(i,:))/norm(RFs_dF(i,:));
end

figure; hold on
plot(TI_sp,TI_dF,'ko')

figure; hold on
plot_scatterRLine(TI_sp,TI_dF)
%% bootstrap test for whisker responsiveness

framesEvoked=int(1):int(2);

for K=1:length(cellPaths)
    load([cellPaths{K},'results.mat'])
    [bootStrapResults,bs_meanDF,medianResponses ] = bootStrap_wrapper_single(...
        results.raw.deltaF,results.byWhisk.traceByStim,10000,sampRateIm,0.5,1,framesEvoked)
    
      pvals=cellfun(@(x)bootStrapResults.(x),whisk,'Uni',1);
    sig_inds_whisk=MultControl(pvals,0.05,'FDR');
    sig_inds(K)=sum(sig_inds_whisk)>0;
end

