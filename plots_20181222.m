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

%% determine which cells are whisker-responsive
sig_inds=zeros(length(cellPaths),1);
for K=1:length(cellPaths)
    pvals=permuteTest_calibration(cellPaths{K},10000,'mean');
    whisk=fieldnames(pvals);
    pvtmp=cellfun(@(x)pvals.(x),whisk);
    pvtrue=MultControl(pvtmp,0.05,'FDR');
    sig_inds(K)=sum(pvtrue)>0;
end
sig_inds=logical(sig_inds);
%%
cellPaths=cellPaths(sig_inds);
%%
colormap gray
for K=1:length(cellPaths)
    %load in reduced data for this cell
    load(strcat(cellPaths{K},'results.mat'));
    
    figure;
%     subplot(1,2,1)
    RF_heatmap(results.byWhisk.traceByStim)
    title(cellPaths{K})    
fnstmp=fieldnames(results.raw.ImMetadata);
       sampRateCA=results.raw.ImMetadata.(fnstmp{1}).sampRateCA;
    sampRateIm=results.raw.ImMetadata.(fnstmp{1}).sampRateIm;

    whisk=results.byWhisk.whisk;   
    
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

end
