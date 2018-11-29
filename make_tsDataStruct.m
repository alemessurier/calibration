function timeSeries=make_tsDataStruct(cellPaths,savePath)

cd(savePath);
    tmp=dir('timeSeries.mat');
    if ~isempty(tmp)
        load(strcat(savePath,'timeSeries.mat'),'timeSeries')
        numApp=size(timeSeries,2);
    else
        numApp=0;
    end
    
for K=1:length(cellPaths)
    load(strcat(cellPaths{K},'imgFileInfo.mat'),'imDir','imFiles');
    load(strcat(cellPaths{K},'results.mat'),'results')
    ROIname=strcat('ROI',num2str(K));
    for J=1:length(imFiles)
        [ImageArray, Metadata] = LoadTIFF_SI5(strcat(imDir,imFiles{J},'.tif'));
        [Metadata] = ReadTIFFHeader_SI5(Metadata);
        ROIMask=define_ROI_masks(ImageArray,{});
        npMask=makeNeuropilMasks(ROIMask);
        
        tmproiTS = get_fluoTimeSeries_OL( ROIMask,ImageArray );
        tmpnpTS= get_fluoTimeSeries_OL( npMask,ImageArray );
        roiTS=tmproiTS.ROI1;
        npTS=tmpnpTS.ROI1;
        
        timeSeries(K+numApp).ROIMask.(imFiles{J})=ROIMask;
        timeSeries(K+numApp).npMask.(imFiles{J})=npMask;
        timeSeries(K+numApp).roiTS.(imFiles{J})=roiTS;
        timeSeries(K+numApp).npTS.(imFiles{J})=npTS;
    end
         timeSeries(K+numApp).path=cellPaths{K};
        save(strcat(savePath,'timeSeries.mat'),'timeSeries')
  
end

    
    
    
    function npMask=makeNeuropilMasks(ROIMask)
      
        g = exp(-(-10:10).^2/2/2^2);
        maskb = conv2(g,g,double(logical(ROIMask)),'same')>.15;                        % dilation for border region around ROIs
        [xi,yi] = meshgrid(1:size(ROIMask,1),1:size(ROIMask,2));
        centroid=regionprops(ROIMask(:,:),'centroid');
            centroids=centroid.Centroid;
            for neuropilrad = 40:5:100
                M = (xi-centroids(1)).^2+(yi-centroids(2)).^2 < neuropilrad^2;    % mask of pixels within the radius
                npMask(:,:) = M.*~maskb;                                          % remove ROIs and border regions
                if nnz(npMask(:,:)) > 4000
                    break
                end
            end
        end
end