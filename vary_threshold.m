function [ false_pos_rate,hit_rate,dPrime ] = vary_threshold( cellPaths )
%VARY_THRESHOLD Wrapper function for THRESHOLD_DETECT - varies thresholds
%used to detect spikes from deltaF/F; calculates dPrime for each value
%
% INPUT
%           'cellPaths'         list of directories containing results.mat file
%                               for each cell
% OUTPUTS
%           'false_pos_rate'    cell array of false positive rates when
%                               threshold and delay period for spike detction are varied
%
%           'hit_rate'          cell array of hit rates when threshold and delay period for spike detction are varied
%
%           'dPrime'            cell array of dPrimes for false_pos_rate
%                               and hit_rate
%
thresh_scale=0:0.01:1; %range of thresholds for deltaF/F
delay_scale=0:0.1:1;  %range of delays from spike times in which to look for threshold crossings

figure; hold on
cmap=morgenstemning(length(thresh_scale)+2);
cmap=cmap(2:end,:);

for j=1:length(delay_scale) % loop through delays
    for i=1:length(thresh_scale) % loop through thresholds
        detect_all=[]; 
        false_detect_all=[];
        for J=1:length(cellPaths) %loop through all recorded cells
            load(strcat(cellPaths{J},'results.mat')); %load in data
            [ detect_all_cp, false_detect_all_cp ] = threshold_detect( results,thresh_scale(i),delay_scale(j));
            detect_all=[detect_all detect_all_cp]; %append results for each cell
            false_detect_all=[false_detect_all false_detect_all_cp];
            clear results
        end
        hit_rate{j}(i)=sum(detect_all)/length(detect_all); %calculate hit rate for each threshold/delay combo
        false_pos_rate{j}(i)=sum(false_detect_all)/length(false_detect_all); %calculate false alarm rate for each threshold/delay combo
        dPrime{j}(i)=dprime1(hit_rate{j}(i),false_pos_rate{j}(i),length(detect_all),length(detect_all)); %calculate d' for each threshold/delay combo
    end
    plot(false_pos_rate{j},hit_rate{j},'Color',cmap(j,:))

end
% figure
end

