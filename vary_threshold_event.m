function [ false_pos_rate,hit_rate,correct_rej_rate, miss_rate,dPrime ] = vary_threshold_event( cellPaths )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
thresh_scale=0:0.1:1;
delay_scale=0:0.1:1;
figure; hold on
cmap=morgenstemning(length(thresh_scale)*5);

for j=1:length(delay_scale)
    for i=1:length(thresh_scale)
        detect_all=[];
        false_detect_all=[];
         misses_all=[];
            correct_rej_all=[];
        for J=1:length(cellPaths)
            load(strcat(cellPaths{J},'results.mat'));
            [ detect_all_cp, false_detect_all_cp,misses_all_cp,correct_rej_all_cp  ] = threshold_detect_event( results,thresh_scale(i),delay_scale(j));
            detect_all=[detect_all detect_all_cp];
            false_detect_all=[false_detect_all false_detect_all_cp];
            misses_all=[misses_all misses_all_cp];
            correct_rej_all=[correct_rej_all correct_rej_all_cp];
            clear results
        end
        hit_rate{j}(i)=sum(detect_all)/(sum(detect_all)+sum(misses_all));
        false_pos_rate{j}(i)=sum(false_detect_all)/(sum(false_detect_all)+sum(correct_rej_all));
        miss_rate{j}(i)=sum(misses_all)/length(misses_all);
        correct_rej_rate{j}(i)=sum(correct_rej_all)/length(correct_rej_all);
        dPrime{j}(i)=dprime1(hit_rate{j}(i),false_pos_rate{j}(i),sum(detect_all)+sum(misses_all),sum(false_detect_all)+sum(correct_rej_all));
    end
    plot(false_pos_rate{j},hit_rate{j},'o','MarkerEdgeColor',cmap(j*4,:))

end
% figure
end

