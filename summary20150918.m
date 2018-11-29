%% Loop through all cells; make templates for Ca transients resulting from 1,2,3,4,and <=5 spikes

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

for K=1:length(cellPaths)
    load(strcat(cellPaths{K},'results.mat'));
    [ template100ms(K) ] = make_spikeCountTemplate( results,0.1 );
    [ template500ms(K) ] = make_spikeCountTemplate( results,0.5 );
    [ template1sec(K) ] = make_spikeCountTemplate( results,1 );
end

sampRateCA=results.raw.ImMetadata.cell3_00001_grch.sampRateCA;
sampRateIm=results.raw.ImMetadata.cell3_00001_grch.sampRateIm;
fnstmp=fieldnames(template100ms);
h=figure; hold on
i=figure; hold on
j=figure; hold on
for J=1:length(fnstmp)
    tmp=arrayfun(@(x)x.(fnstmp{J}),template100ms,'Uni',0);
    spTemp100.(fnstmp{J})=vertcat(tmp{:});
    meanTmp100=mean(spTemp100.(fnstmp{J}),1);
    figure(h);
    plot((1:length(meanTmp100))/sampRateIm -0.5,meanTmp100);

    tmp=arrayfun(@(x)x.(fnstmp{J}),template500ms,'Uni',0);
    spTemp500.(fnstmp{J})=vertcat(tmp{:});
    meanTmp500=mean(spTemp500.(fnstmp{J}),1);
    figure(i);
    plot((1:length(meanTmp500))/sampRateIm -0.5,meanTmp500);

    tmp=arrayfun(@(x)x.(fnstmp{J}),template1sec,'Uni',0);
    spTemp1sec.(fnstmp{J})=vertcat(tmp{:});
    meanTmp1sec=mean(spTemp1sec.(fnstmp{J}),1);
    figure(j);
    plot((1:length(meanTmp1sec))/sampRateIm -0.5,meanTmp1sec);
end

title('Mean response to N spikes/0.1 sec')
legend({'N=1','N=2','N=3','N=4','N>=5'})

[maxDF100,maxInds100]=cellfun(@(x)max(mean(spTemp100.(x),1)),fnstmp,'Uni',1)
[maxDF500,maxInds500]=cellfun(@(x)max(mean(spTemp500.(x),1)),fnstmp,'Uni',1)
[maxDF1sec,maxInds1sec]=cellfun(@(x)max(mean(spTemp1sec.(x),1)),fnstmp,'Uni',1)
figure
plot(maxInds100)
hold on
plot(maxInds500)
plot(maxInds1sec)

figure
plot(maxDF100)
hold on
plot(maxDF500)
plot(maxDF1sec)
%%
 cellPaths={'J:\reduced\20150817\cell2\',...
    'J:\reduced\20150819\cell2\',...
    'J:\reduced\20150824\cell1\',...
    'J:\reduced\20150824\cell2\',...
    'J:\reduced\20150824\cell4\',...
    'J:\reduced\20150825\cell2\',...
    'J:\reduced\20150825\cell1\',...
    'J:\reduced\20150827\cell1\',...
    'J:\reduced\20150827\cell2\',...
    'J:\reduced\20150827\cell3\'...
    };

%%
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

all_templates=[];
all_inds.cumSumSp=[];
all_inds.ISI_Pre=[];
all_inds.ISI_Post=[];
all_bsTemplates=[];

for K=1:length(cellPaths)
    load(strcat(cellPaths{K},'results.mat'));
    [ templates{K},templateInds(K),bootStrap_templates{K}] = make_spikeCountTemplateV3( results );
    all_bsTemplates=[all_bsTemplates; bootStrap_templates{K}];
    all_templates=[all_templates; templates{K}];
    all_inds.cumSumSp=[all_inds.cumSumSp; templateInds(K).cumSumSp];
    all_inds.ISI_Pre=[all_inds.ISI_Pre; templateInds(K).ISI_pre];
    all_inds.ISI_Post=[all_inds.ISI_Post; templateInds(K).ISI_post];

end

%%
firstSpInds=all_inds.ISI_Pre>2;
templates_firstSp=templates_bs(firstSpInds);
templates_firstSp_num=horzcat(templates_firstSp{:})';
templates_bs=cellfun(@(x)x-mean(x(1:3)),all_templates,'Uni',0);
cumSum_firstSp=all_inds.cumSumSp(firstSpInds,:);
[maxDFs,maxDFinds]=cellfun(@max,templates_firstSp);
maxDFs_bs=max(all_bsTemplates,[],2);
maxDFs_bs=sort(maxDFs_bs,'ascend');
bs_percentile_95=maxDFs_bs(floor(0.95*length(maxDFs_bs)))
cmap=morgenstemning(25);
figure; hold on
for i=1:20
    spikeNums=cumSum_firstSp(:,i*1000);
%     uniqueNums=unique(spikeNums);
    meanMaxDF=zeros(1,10);%length(uniqueNums));
    for j=1:10;%length(uniqueNums)
%         num=uniqueNums(j);
        meanMaxDF(j)=mean(maxDFs(spikeNums==j));
    end
    plot(1:10,meanMaxDF,'Color',cmap(i,:))
    legend_labels{i}=strcat('t=',num2str(i/10),'sec');
end
hline(bs_percentile_95);
legend(legend_labels);
xlabel('num spike in t sec')
ylabel('mean dF/F peak value')

%%
firstSpMeans=mean(templates_firstSp(:,4:10),2);
figure
imagesc(all_inds.cumSumSp(firstSpInds,:))
figure
imagesc(templates_firstSp)
tmp=all_inds.cumSumSp(firstSpInds,:);

figure
plot(firstSpMeans,all_inds.cumSumSp(firstSpInds,10000),'ko')


cumSum_firstSp=all_inds.cumSumSp(firstSpInds,:);
[numSpikesSort,numSpikesInds]=sort(cumSum_firstSp(:,500),'descend');
templates_sortNumSpikes=templates_firstSp(numSpikesInds);
templates_sortNumSpikes=horzcat(templates_sortNumSpikes{:})';
figure
imagesc(templates_sortNumSpikes)
figure
surf(templates_sortNumSpikes)
%%

% 



templates_bs=cellfun(@(x)x-mean(x(1:3)),all_templates,'Uni',0);
template_means=cellfun(@(x)mean(x(4:10)),templates_bs,'Uni',1);
[means,IndsMean]=sort(template_means,'descend');
templatesSortMean=templates_bs(IndsMean);
templatesSortMean=horzcat(templatesSortMean{:})';
figure; imagesc(templatesSortMean)

template_max=cellfun(@(x)max(x(4:10)),templates_bs,'Uni',1);
[maxDF,IndsMax]=sort(template_max,'descend');
templatesSortMax=templates_bs(IndsMax);
templatesSortMax=horzcat(templatesSortMax{:})';
figure; imagesc(templatesSortMax)
surf(templatesSortMax)
% make heat map of transients sorted by ISI before spike

[~,ISI_B]=sort(all_inds.ISI_B,'descend');

templates_ISI_B=templates_bs(ISI_B);
templates_ISI_B=horzcat(templates_ISI_B{:})';


% make heat map of transients sorted by ISI after spike

[~,ISI_A]=sort(all_inds.ISI_A,'descend');

templates_ISI_A=templates_bs(ISI_A);
templates_ISI_A=horzcat(templates_ISI_A{:})';
figure; imagesc(templates_ISI_A)

%%

firstSpInds=all_inds.ISI_B>2;
templates_firstSp=templates_bs(firstSpInds);
[firstSpISIA,ISI_A]=sort(all_inds.ISI_A(firstSpInds),'descend');
templates_ISI_A=templates_firstSp(ISI_A);
templates_ISI_A=horzcat(templates_ISI_A{:})';
figure; imagesc(templates_ISI_A)
meanDFs_ISI_A=mean(templates_ISI_A(:,4:10),2);
figure
plot(firstSpISIA,meanDFs_ISI_A)

template_means=cellfun(@(x)mean(x(4:end)),templates_firstSp,'Uni',1);
[means,IndsMean]=sort(template_means,'descend');
templatesSortMean=templates_firstSp(IndsMean);
templatesSortMean=horzcat(templatesSortMean{:})';
figure; surf(templatesSortMean)

template_max=cellfun(@(x)max(x(4:end)),templates_firstSp,'Uni',1);
[maxDF,IndsMax]=sort(template_max,'descend');
templatesSortMax=templates_firstSp(IndsMax);
templatesSortMax=horzcat(templatesSortMax{:})';
figure; 
surf(templatesSortMax)


%% 

%%
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

all_templates=[];
all_inds.ISI_Post=[];
all_inds.ISI_Pre=[];

for K=1:length(cellPaths)
    load(strcat(cellPaths{K},'results.mat'));
    [ templates,templateInds] = make_spikeCountTemplateV2( results );
    
        all_templates=[all_templates; templates];
    all_inds.cumSumSp=[all_inds.cumSumSp; templateInds.cumSumSp];
    all_inds.ISI_Pre=[all_inds.ISI_Pre; templateInds.ISI_pre];

end

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

int(1)=ceil(0.5*sampRateIm);
int(2)=ceil(0.5*sampRateIm)+round(0.5*sampRateIm);
colormap gray
areaSp=zeros(1,9);
areaDF=zeros(1,9);
for K=1:length(cellPaths)
    load(strcat(cellPaths{K},'results.mat'));
    whisk=results.byWhisk.whisk;
    figure
    subplot(1,2,1)
    responseVecSp=mean(results.spikeTuning.pSpike,2);
    [~,areaSp(K)]=calculateCentroid(responseVecSp,0.5,1);
    title('spikes')
    responseVec_dF=cellfun(@(x)mean(mean(results.byWhisk.traceByStim.(x)(:,int(1):int(2)))),whisk,'Uni',1);
    axis square
    subplot(1,2,2)
    [~,areaDF(K)]=calculateCentroid(responseVec_dF,0.5,1);
   axis square
    colormap gray
       title('deltaF')
    
    RF_proj(K)=acos(dot(responseVecSp,responseVec_dF)/(norm(responseVecSp)*norm(responseVec_dF)))*(180/pi);

end

tmp=cellPaths(areaSp>areaDF)

