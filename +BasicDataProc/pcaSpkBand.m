function [pcaProc] = pcaSpkBand(spkSmoothed, bandPower, lbls)
%function to fit in procDataX(CO, GT, HeadM, etc)
%   gets a pca and analyzes the peaks

%Inputs
% spkSmoothed- spike smoothed data using Analysis.BasicDataProc.spikeRate output
% bandPower- classic or other band, smoothed power analysis using Analysis.BasicDataProc.dataPrep output
% lbls- typically can input fields(bandPower)



for nn=1:length(lbls)
    for ii=1:size(bandPower,2)
        [~, spikeTemp(:, :, ii)]=pca(spkSmoothed{ii}); %time x pca/ch x trial
        if nn==1
            pcaProc(ii).spikeSmooth=spikeTemp(:, :, ii);
            pcaProc(ii).(strcat('spikeSmooth', 'Peaks'))=Analysis.BasicDataProc.specChange(spikeTemp(:,:,ii)); %find the peaks
        end
        [~, bandTemp(:,:,ii)]=pca(squeeze(bandPower{ii}(:,nn,:))); %time x pca/ch x trial
        pcaProc(ii).(lbls{nn})=bandTemp(:,:,ii);
        pcaProc(ii).(strcat(lbls{nn}, 'Peaks'))=Analysis.BasicDataProc.specChange(bandTemp(:,:,ii)); %find the peaks
    end
    
    pcaProc(1).(strcat(lbls{nn}, 'TrialMean'))=nanmean(bandTemp,3); %mean across all trials
    pcaProc(1).(strcat(lbls{nn}, 'SE'))=nanstd(bandTemp, [], 3)/sqrt(size(bandTemp,3));
    if nn==length(lbls)
        pcaProc(1).(strcat('spikeSmooth', 'TrialMean'))=nanmean(spikeTemp,3);
        pcaProc(1).(strcat('spikeSmooth', 'TrialSE'))=nanstd(spikeTemp, [], 3)/sqrt(size(spikeTemp,3));
    end    
end


end

