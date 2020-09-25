function [dataHoldMove] = procDataCO(R, ns5File, varargin)
%procDataGridT Process the BCI data to break it up into times when it is
%hovering and when it isnt for the center out (GT) task.
%  R details - each row is a ms
% Input
%     R - R struct file, 
%     ns5File - the name of the file path for the file you want
%     
% Output
%     dataHoldMove - has 6 outputs
%%TIMES
% Remember: currently epoch is 300ms
% Hold=300ms before hover period is achieved (11)
% Move=300ms starting 200 ms after the target appears
% TransHold=300ms starting 200ms prior to hold state (hold is 500 ms, so it's 700ms before click state of 11) into 100 ms of hold state
% TransMove=300ms starting 100ms before the target appears (state 3)
%%
% TO DO:
%    
% 
%     Example:
% 
% %Task open
% ns5file='C:\Users\dankr\Documents\Data\121218_b7_NPTL_COData\LateralGrid_121218_b7\7_cursorTask_Complete_t5_bld(007)008.ns5';
% load('C:\Users\dankr\Documents\Data\121218_b7_NPTL_COData\cursorTaskR_T5.mat')
% R=R(1:64);
%         
%         [dataHoldMove]=Analysis.NPTL.procDataGridT(R, ns5file, 'Spectrogram', true);
%
% Blackrock gridlayout
%   NaN     2     1     3     4     6     8    10    14   NaN
%     65    66    33    34     7     9    11    12    16    18
%     67    68    35    36     5    17    13    23    20    22
%     69    70    37    38    48    15    19    25    27    24
%     71    72    39    40    42    50    54    21    29    26
%     73    74    41    43    44    46    52    62    31    28
%     75    76    45    47    51    56    58    60    64    30
%     77    78    82    49    53    55    57    59    61    32
%     79    80    84    86    87    89    91    94    63    95
%    NaN    81    83    85    88    90    92    93    96   NaN

% THINGS TO DO FOR PROCESSING: FIRST PASS. SO MAIN THING TO DO IS TO PLOT THEM ALL
%     RUN A CORRELATION AT THE DIFFERENT ELECTRODES ACROSS THE GRID
%         MAKE SURE TO OUTPUT THAT RAW DATA TO RUN IT BETWEEN THE GRIDS AS WELL
%     LOOK AT THE CHANGES IN EACH BAND DURING CHANGES IN FIRING RATE.
[varargin, data] = util.argkeyval('data',varargin, []); %check if the data has already been extracted    
[varargin, fs] = util.argkeyval('fs',varargin, 30000); %original sample rate
[varargin, downSample] = util.argkeyval('downSample',varargin, 2000); %the new downsample rate

[varargin, epoch] = util.argkeyval('epoch',varargin, 500); %the epoch for all times, in ms
[varargin, spkWin] = util.argkeyval('spkWin',varargin, 20); %the window to average spikes over, in ms
%These are the times to take around the events (state==11, click state,
%after holding for 500ms, and state==3, appearance of a new target). These
%are set up so for instance transMovePre starts 100ms prior to target
%appearing, under the assumption that shows the transition to movement.
[varargin, holdPre] = util.argkeyval('holdPre',varargin, 0); %time before transition to 11, a state of click, hovering for 500ms, in ms
[varargin, transHoldPre] = util.argkeyval('transHoldPre',varargin, -400); %time before transition to 11, a state of click, hovering for 500ms, in ms
[varargin, movePre] = util.argkeyval('movePre',varargin, 0); %time after a target appears, transition to state 3
[varargin, transMovePre] = util.argkeyval('transMovePre',varargin, -100); %time before a target appears, transition to state 3

[varargin, bufferT] = util.argkeyval('bufferT',varargin, 300); %the buffer window in ms to process data around your desired time of interest (so will take bufferT s before and bufferT after the epoch). This is to allow sufficient ramp up and plotting before and after
[varargin, channel] = util.argkeyval('channel',varargin, [1 96]); %which channels you want, default is whole grid
[varargin, trialNum] = util.argkeyval('channel',varargin, [1, size(R,2)]); %which trials you want, default is all.
[varargin, perc] = util.argkeyval('channel',varargin, 95); %what percentile of the most active channels you want. THIS LIMITS HOW MANY CHANNELS ARE INCLUDED, HIGHER PERCENTILE MEANS LESS CHANNELS ARE DISPLAYED
[varargin, grid] = util.argkeyval('grid', varargin, 1); %1 = latera, 0 = medial.

[varargin, Spectrogram] = util.argkeyval('Spectrogram',varargin, true); %if doing Spectrogram or not

[varargin, plotData] = util.argkeyval('plotData',varargin, true); %amount to plot before the epoch.  Will error if longer than the bufferT.  Only applies to the spectrogram and bandpassed dataolmk,
[varargin, plotPre] = util.argkeyval('plotPre',varargin, 100); %amount to plot before the epoch.  Will error if longer than the bufferT.  Only applies to the spectrogram and bandpassed dataolmk,
[varargin, plotPost] = util.argkeyval('plotPost',varargin, 400); %amount to plot after the epoch.

%% PREP ACTIVITY
%convert R to a struct PROBABLY THIS TO UNWIELDY AND WILL NEED TO BREAK
%R UP
if ~isa(R, 'struct')
    R=[R{:}];
end

%% OPEN DATA
%convert to readable for the channels
ch=['c:0', num2str(channel(1,1)), ':', num2str(channel(1,2))];
if isempty(data)
    ns5=openNSx('read', ns5File, ch); %Open the file
    data=double(ns5.Data)';
end

%make sure data is columns=channels and rows are time
if size(data, 1)<size(data,2)
    data=data';
end

%gridlayout per idx to channel FOR NOW, ONLY ADJUSTING THE CHANNELS IN
%CORRELATION, EASIER TO KEEP TRACK OF WHICH CHANNEL IS WHICH
idxch2el=[78,88,68,58,56,48,57,38,47,28,37,27,36,18,45,17,46,8,35,16,24,7,26,6,25,5,15,4,14,3,13,2,77,67,76,66,75,65,74,64,73,54,63,53,72,43,62,55,61,44,52,33,51,34,41,42,31,32,21,22,11,23,10,12,96,87,95,86,94,85,93,84,92,83,91,82,90,81,89,80,79,71,69,70,59,60,50,49,40,39,30,29,19,20,1,9];


lblN={'DeltaOrSpike', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High_Gamma'};
lblP={'MaxStart', 'MaxStartCILow', 'MaxStartCIHigh', 'MaxPeak', 'MaxPeakCILow',...
    'MaxPeakCIHigh', 'MinStart', 'MinStartCILow', 'MinStartCIHigh', 'MinPeak',...
    'MinPeakCILow',    'MinPeakCIHigh'};
idxL=1;
for ii=1:length(lblN)
    for jj=1:length(lblP)
        Peak(idxL).Runs___________=strcat(lblN{ii}, lblP{jj});
        idxL=idxL+1;
    end
end


%create an iti, idle time is 3 seconds prior to start, end .5 seconds early
%to prevent anticipation
%THE TRIAL IS NOT ALWAYS THE FIRST ONE IF IT'S .isSuccessful ONLY, THIS IS
%A WORK AROUND
meanTrialTime=floor(mean(double([R.trialLength])));
firstTrial=double(R(1).trialNum);
preTrialTime=meanTrialTime*firstTrial*30; %subtract this to get to around the beginning, add another second to be sure, multiply by 30 fs
itiD=data(double(R(1).firstCerebusTime(1,1))-(30000*3)-300-preTrialTime:double(R(1).firstCerebusTime(1,1))-300-preTrialTime, :);%also subtract the 10ms offset, or 300 samples


dataClassBand=[];
[specTempHoldI, ~,~,~, dataClassBand] =Analysis.BasicDataProc.dataPrep(Analysis.BasicDataProc.downSample(itiD, fs, downSample), 'Spectrogram', Spectrogram, 'doBandFilterBroad', true, 'dataClassBand', dataClassBand, 'itiProc', true);
if isempty(dataClassBand)
    for ii=3:length(lblN)
        dataClassBand.(lblN{ii})=[];
    end
end

lblB=fields(dataClassBand);%for naming things later

tI=specTempHoldI.dataSpec.tplot;
fI=specTempHoldI.dataSpec.f;


%Spike times to set up for plotting
tSpkHold=linspace(holdPre-epoch, 0, floor(abs(epoch)/spkWin)); %t=0 is state 11 starting, prior to click to click (state 3 to 11, 11 being click state )
tSpkTransHold=linspace(transHoldPre-epoch, transHoldPre, floor(abs(epoch)/spkWin)); %t=0 is hold state
tSpkMove=linspace(movePre, movePre+epoch, floor(abs(epoch)/spkWin)); %t=0 is when the new target appears
tSpkTransMove=linspace(transMovePre, transMovePre+epoch, floor(abs(epoch)/spkWin)); %t=0 is when a new target appears
tGen=linspace(0, epoch, floor(abs(epoch)/spkWin)); %generic 300ms for easy comparison of epochs

%% gather the data up
holdData=[]; moveData=[];

idx1=2; trl=1; fl=1;
rState=double([R.state]);
rCP=double([R.cursorPosition]); %1 and 3 are position, 2 and 4 are velocity for grid task, and it's derivative of cp in grid task or row 1:2 of xp
rFCT=double([R.firstCerebusTime])+10*30; %There is a 10ms offset from ns5 and firstCerebusTime, unknown why, but verified. difference between row 1 and 2 is 6ms, lastcerebrustime is 29 and 35ms respectively for row 1/2
rMinA=double([R.minAcausSpikeBand]);
if isempty(cell2mat(strfind(fields(R), 'spikeRaster')))
    threshMultiplier = -4.5;
    RMSvals = channelRMS( R );
    thresholds = RMSvals.*threshMultiplier;
    thresholds = thresholds'; % make column vector
    for iTrial = 1 : numel( R )
        R(iTrial).spikeRaster = R(iTrial).minAcausSpikeBand < repmat( thresholds, 1, size( R(iTrial).minAcausSpikeBand, 2 ) );
    end
end
rSpkAll=double([R.spikeRaster]);
%choose your grid to evaluate, depends on the ns5. 
if grid
    rSpk=rSpkAll(1:96,:);
    rMinA=rMinA(1:96,:);
else
    rSpk=rSpkAll(97:192, :);
    rMinA=rMinA(97:192,:);
end



%% spikes, all channels
%smooth the spikes
[spk, spkRate, spkRateSmooth, tm] = Analysis.BasicDataProc.spikeRate(rSpk, 'win', spkWin);
% Outputs:l,
%         spk-  converts raster to NaNs and 1s for plotting
%         spkRate- spikes summed over a window
%         spkRateSmooth- smoothed spikes with convolve
%         tm- a time vector for plotting to match spikes
spikeCh.spk=spk; %flip them all to match the lfp data, so columns are channels
spikeCh.spkRate=spkRate;
spikeCh.spkRateSmooth=spkRateSmooth;
spikeCh.tm=tm;


%get the states hovering
idxH=1; idxM=1; idxCheck=1; idxT=1;
for ii=2:length(rState)-1 %first trial is different than the rest
    if rState(ii)==5 && rState(ii-1)==4 %5 is the click state, 4 is the hover state
        rHoldEnd(idxH,1)=ii;
        holdCP{idxH}=rCP(:,ii+holdPre-epoch:ii);
        holdSpk{idxH}=spk(ii+holdPre-epoch:ii, :); %take all the channels at once
        holdFCT(idxH, :)=rFCT(1,ii); 
        holdData{idxH}=data(holdFCT(idxH, 1)+((holdPre-epoch-bufferT)*30):holdFCT(idxH, end)+((bufferT)*30),:); %add the buffer for processing and plotting
        tempSt=find(tm>=ii+holdPre-epoch); tempEnd=find(tm<=ii);
        holdSmooth{idxH}=spkRateSmooth(tempSt(1):tempEnd(end),:); 
        %get spikes averaged across the channels, for each trial
        tempSt=[]; tempEnd=[];
        
        idxH=idxH+1;
    elseif rState(ii)==2 && rState(ii-1)~=2  %2 is when the target first appears,
        rMoveSt(idxM, 1)=ii;
        moveCP{idxM}=rCP(:,ii+movePre:ii+movePre+epoch); %after movepre ms after the target appears to get into a move state
        moveSpk{idxM}=spk(ii+movePre:ii+movePre+epoch, :);
        moveFCT(idxM, :)=rFCT(1,ii);
        moveData{idxM}=data(moveFCT(idxM, 1)+((movePre-bufferT)*30):moveFCT(idxM, end)+((movePre+epoch+bufferT)*30), :); %add the buffers for processing
        tempSt=find(tm>=ii+200); tempEnd=find(tm<=ii+200+epoch);
        moveSmooth{idxM}=spkRateSmooth(tempSt(1):tempEnd(end), :);        
        tempSt=[]; tempEnd=[];
        

        idxM=idxM+1; 
    elseif rState(ii+1)-rState(ii)==9 || rState(ii+1)-rState(ii)<=-8 %check for any weird trials
        rCheck(1,idxCheck)=ii;
    end
end


% find the channels most active, then choose between highest spike rates
% (mean) or highest variability (std)
if size(holdSmooth,2)==96 || size(holdSmooth,2)==96*2 %check which matrix position channels are in, should be time x channels and position should be 2
    chPosition=2;
else
    chPosition=1;
end
[msH] = Analysis.BasicDataProc.activeSpkCh(holdSmooth, 'perc', perc, 'chPosition', chPosition);
[msM] = Analysis.BasicDataProc.activeSpkCh(moveSmooth, 'perc', perc, 'chPosition', chPosition);


%which one to pick to look at the channels
spikeCh.meanstdSpkPerc=msM;
%choose the modDSPk, can also choose mean and ste changes, but likely will
%pick up noisey channels
ChSpk=msM.modDSpk;

%% hold
%average across all channels, for each trial
TT=tic;
for ii=1:size(holdData,2)
    tt=tic;
    [holdDataAllTemp, params, holdCBAllTemp]=Analysis.BasicDataProc.dataPrep(Analysis.BasicDataProc.downSample(holdData{ii}, fs, downSample), 'Spectrogram', Spectrogram, 'dataClassBand', dataClassBand, 'doBandFilterBroad', true);
    toc(tt)
    holdCBAll{ii}=holdCBAllTemp(downSample*bufferT/1000+1:end-downSample*bufferT/1000,:,:); %save the band outputs   
    holdDataAllCh(:,:,ii)=nanmean(holdDataAllTemp.dataSpec.data,3); %mean across all channels    
    holdBandDataAllCh(:,:,ii)=nanmean(holdCBAllTemp,3); %mean across all channels
    holdSmoothAllCh(:,ii)=nanmean(holdSmooth{ii},2);   %mean across all channels 
    %%find maxima and minima of the spikes and band data for each trial,
    %%this is too compare within each trial.
    holdBandPeaksAllCh(:,:,ii)=Analysis.BasicDataProc.specChange(holdBandDataAllCh(downSample*bufferT/1000+1:end-downSample*bufferT/1000,:,ii)); %remove the buffered ends for this analysis.
    holdSpkPeaksAllCh(:,ii)=Analysis.BasicDataProc.specChange(holdSmoothAllCh(:,ii));
    for jj=1:length(ChSpk)
        bychHoldSpk{jj}(:, ii)=holdSmooth{ii}(:, ChSpk(jj)); %cells by channel, with trials by spike averages        
        bychHoldData{jj}(:,:,ii)=holdDataAllTemp.dataSpec.data(:,:,ChSpk(jj)); %load a 3d matrix for time x freq x trial by channel for spectral data
        bychHoldBandData{jj}(:,:,ii)=holdCBAllTemp(:,:,ChSpk(jj));        
    end
end

%% move
toc(TT);
for ii=1:size(moveData,2)
    tt=tic;
    [moveDataAllTemp, params, moveCBAllTemp]=Analysis.BasicDataProc.dataPrep(Analysis.BasicDataProc.downSample(moveData{ii}, fs, downSample), 'Spectrogram', Spectrogram, 'dataClassBand', dataClassBand, 'doBandFilterBroad', true);
    toc(tt)
    moveCBAll{ii}=moveCBAllTemp(downSample*bufferT/1000+1:end-downSample*bufferT/1000,:,:);    
    moveDataAllCh(:,:,ii)=nanmean(moveDataAllTemp.dataSpec.data,3); %mean across all channels    
    moveBandDataAllCh(:,:,ii)=nanmean(moveCBAllTemp,3); %mean across all channels
    moveSmoothAllCh(:,ii)=nanmean(moveSmooth{ii},2);   %mean across all channels 
    %%find maxima and minima of the spikes and band data for each trial
    moveBandPeaksAllCh(:,:,ii)=Analysis.BasicDataProc.specChange(moveBandDataAllCh(downSample*bufferT/1000+1:end-downSample*bufferT/1000,:,ii)); %remove the buffered ends for this analysis.
    moveSpkPeaksAllCh(:,ii)=Analysis.BasicDataProc.specChange(moveSmoothAllCh(:,ii));
    for jj=1:length(ChSpk)
        bychMoveSpk{jj}(:, ii)=moveSmooth{ii}(:, ChSpk(jj)); %cells by channel, with trials by spike averages        
        bychMoveData{jj}(:,:,ii)=moveDataAllTemp.dataSpec.data(:,:,ChSpk(jj)); %load a 3d matrix for time x freq x trial by channel for spectral data
        bychMoveBandData{jj}(:,:,ii)=moveCBAllTemp(:,:,ChSpk(jj));        
    end
end


%% PCA
for nn=1:length(lblB)
    for ii=1:size(moveCBAll,2)
        [~, spikeTemp(:, :, ii)]=pca(moveSmooth{ii}); %time x pca/ch x trial
        if nn==1
            movePCA(ii).spikeSmooth=spikeTemp(:, :, ii);
            movePCA(ii).(strcat('spikeSmooth', 'Peaks'))=Analysis.BasicDataProc.specChange(spikeTemp(:,:,ii)); %find the peaks
        end
        [~, bandTemp(:,:,ii)]=pca(squeeze(moveCBAll{ii}(:,nn,:))); %time x pca/ch x trial
        movePCA(ii).(lblB{nn})=bandTemp(:,:,ii);
        movePCA(ii).(strcat(lblB{nn}, 'Peaks'))=Analysis.BasicDataProc.specChange(bandTemp(:,:,ii)); %find the peaks
    end
    
    movePCA(1).(strcat(lblB{nn}, 'TrialMean'))=nanmean(bandTemp,3); %mean across all trials
    movePCA(1).(strcat(lblB{nn}, 'SE'))=nanstd(bandTemp, [], 3)/sqrt(size(bandTemp,3));
    if nn==length(lblB)
        movePCA(1).(strcat('spikeSmooth', 'TrialMean'))=nanmean(spikeTemp,3);
        movePCA(1).(strcat('spikeSmooth', 'TrialSE'))=nanstd(spikeTemp, [], 3)/sqrt(size(spikeTemp,3));
    end    
end
for nn=1:length(lblB)
    for ii=1:size(holdCBAll,2)
        [~, spikeTemp(:, :, ii)]=pca(holdSmooth{ii}); %time x pca/ch x trial
        if nn==1
            holdPCA(ii).spikeSmooth=spikeTemp(:, :, ii);
            holdPCA(ii).(strcat('spikeSmooth', 'Peaks'))=Analysis.BasicDataProc.specChange(spikeTemp(:,:,ii)); %find the peaks
        end
        [~, bandTemp(:,:,ii)]=pca(squeeze(holdCBAll{ii}(:,nn,:))); %time x pca/ch x trial
        holdPCA(ii).(lblB{nn})=bandTemp(:,:,ii);
        holdPCA(ii).(strcat(lblB{nn}, 'Peaks'))=Analysis.BasicDataProc.specChange(bandTemp(:,:,ii)); %find the peaks
    end
    
    holdPCA(1).(strcat(lblB{nn}, 'TrialMean'))=nanmean(bandTemp,3); %mean across all trials
    holdPCA(1).(strcat(lblB{nn}, 'SE'))=nanstd(bandTemp, [], 3)/sqrt(size(bandTemp,3));
    if nn==length(lblB)
        holdPCA(1).(strcat('spikeSmooth', 'TrialMean'))=nanmean(spikeTemp,3);
        holdPCA(1).(strcat('spikeSmooth', 'TrialSE'))=nanstd(spikeTemp, [], 3)/sqrt(size(spikeTemp,3));
    end    
end

%%
tmBF=1000/downSample; %size of each bin
tmSF=epoch/size(holdSmooth{1},1); %size of each bin
bandNum=size(bychMoveBandData{1},2);
%go through each channel and pull the spike/lfp data
for jj=1:length(ChSpk) %run through relevant channels    
    chstr=['ch', num2str(ChSpk(jj))];
    perCh.(chstr).meansteSpkMove(:,1)=nanmean(bychMoveSpk{jj},2);   
    perCh.(chstr).meansteSpkMove(:,2)=nanstd(bychMoveSpk{jj}, [], 2)/sqrt(size(bychMoveSpk{1},2));
    peaksTemp=Analysis.BasicDataProc.specChange(bychMoveSpk{jj}); %get the peaks   
    ciTemp=Analysis.BasicDataProc.maxminCI(peaksTemp, 'binSize', tmSF); %get the means and CI, remember this is mean and CI of the peaks themselves, not a mean/CI of the signal and then finding the peak of that. 
    for mm=1:length(ciTemp); Peak(mm).(strcat(chstr, 'Move', 'Spike'))=ciTemp(mm); end %load into peaks
    
    perCh.(chstr).meansteBandMove(:, 1:bandNum)=nanmean(bychMoveBandData{jj},3);
    perCh.(chstr).meansteBandMove(:, bandNum+1:bandNum*2)=nanstd(bychMoveBandData{jj}, [] ,3)/sqrt(size(bychMoveBandData{jj},3));    
    peaksTemp=Analysis.BasicDataProc.specChange(bychMoveBandData{jj}(downSample*bufferT/1000+1:end-downSample*bufferT/1000,:,:)); %remove the buffered ends for this analysis.
    ciTemp=Analysis.BasicDataProc.maxminCI(peaksTemp, 'binSize', tmBF);
    for  ii=1:size(ciTemp); Peak(ii).(strcat(chstr, 'Move', 'Band'))=ciTemp(ii);    end
    
    perCh.(chstr).meansteSpkHold(:,1)=nanmean(bychHoldSpk{jj}, 2);
    perCh.(chstr).meansteSpkHold(:,2)=nanstd(bychHoldSpk{jj}, [], 2)/sqrt(size(bychHoldSpk{1},2));
    peaksTemp=Analysis.BasicDataProc.specChange(bychHoldSpk{jj});    
    ciTemp=Analysis.BasicDataProc.maxminCI(peaksTemp,  'binSize', tmSF); %get the CI 
    for mm=1:length(ciTemp); Peak(mm).(strcat(chstr, 'Hold', 'Spike'))=ciTemp(mm); end %load into peaks
        
    perCh.(chstr).meansteBandHold(:, 1:bandNum)=nanmean(bychHoldBandData{jj},3);
    perCh.(chstr).meansteBandHold(:, bandNum+1:bandNum*2)=nanstd(bychHoldBandData{jj},[] ,3)/sqrt(size(bychHoldBandData{jj},3));
    peaksTemp=Analysis.BasicDataProc.specChange(bychHoldBandData{jj}(downSample*bufferT/1000+1:end-downSample*bufferT/1000,:,:)); %remove the buffered ends for this analysis.
    ciTemp=Analysis.BasicDataProc.maxminCI(peaksTemp, 'binSize', tmBF);
    for  ii=1:size(ciTemp); Peak(ii).(strcat(chstr, 'Hold', 'Band'))=ciTemp(ii);    end
    
    perCh.(chstr).specMove=nanmean(bychMoveData{jj}, 3);
    perCh.(chstr).specHold=nanmean(bychHoldData{jj}, 3);
    
    %when ready to run the clusters
    %[ mnd1, mnd2, ~, ~, sigclust, rlab ] = stats.cluster_permutation_Ttest_gpu3d( S, Sbase, 'alph', 0.00024 );

end

moveDataAllTrAllCh=nanmean(moveDataAllCh,3); %across all channels and trials
moveBandDataAllTrAllCh=nanmean(moveBandDataAllCh, 3); %across all channels and all trials
%find the peaks and CI
peaksTemp=Analysis.BasicDataProc.maxminCI(moveBandPeaksAllCh, 'binSize', tmBF);
for  ii=1:size(peaksTemp)
    Peak(ii).moveBandAllChAllTr=peaksTemp(ii);
end
peaksTemp=Analysis.BasicDataProc.maxminCI(moveSpkPeaksAllCh, 'binSize', tmSF); %across all channels and all trials
for  ii=1:size(peaksTemp)
    Peak(ii).moveSpkAllChAllTr=peaksTemp(ii);
end


holdDataAllTrAllCh=nanmean(holdDataAllCh,3); %across all channels and trials
holdBandDataAllTrAllCh=nanmean(holdBandDataAllCh, 3); %across all channels and all trials
%find the peaks and CI
peaksTemp=Analysis.BasicDataProc.maxminCI(holdBandPeaksAllCh, 'binSize', tmBF);
for  ii=1:size(peaksTemp)
    Peak(ii).holdBandAllChAllTr=peaksTemp(ii);
end
holdSmoothAllTrAllCh=Analysis.BasicDataProc.maxminCI(holdSpkPeaksAllCh, 'binSize', tmSF); %across all channels and all trials
for  ii=1:size(holdSmoothAllTrAllCh)
    Peak(ii).holdSpkAllChAllTr=holdSmoothAllTrAllCh(ii);
end



if plotData
    N=100;
    C=linspecer(N);
    Analysis.NPTL.PLOTpeakCI; %script to run the 
    
   
        if Spectrogram;  win=params.win(1)*1000; %get the window of the spectrogram
        else; win=0; end;
        fDiff=holdDataAllTemp.dataSpec.f;
        tpDiffFull=(holdDataAllTemp.dataSpec.tplot)*1000-bufferT; %convert to ms, adjust 0 with the buffer, and finally move the spectrogram timing to the beginning of the window rather than the middle.
        specPlotIdx=tpDiffFull>-100 & tpDiffFull<400; %only plot -100 to 400 for an epoch of 300. There are large shifts at the beginning and end (although ramp ups should be dealt with in the data prep)
        tpDiffforPlot=tpDiffFull(specPlotIdx);
        tpHold=tpDiffforPlot-epoch-holdPre; %300ms before hover period ended 
        tpMove=tpDiffforPlot+movePre; %200ms after the target appears
        
        %% plotting
        
       
        colorSh=[2, 15, 62, 98];
        colorShBand=[2, 12, 41, 96];
        
        tSpk=linspace(0, epoch, epoch/spkWin);
        tspkHold=tSpk-epoch-holdPre; %300ms before hover period ended 
        tspkMove=tSpk+movePre; %200ms after the target appears
       
        tBand=linspace(0, epoch+bufferT*2, size(bychHoldBandData{1},1))-bufferT; %include the buffers, than adjust 0 to the buffer
        bandPlotIdx=tBand>-plotPre & tBand<plotPost; %only plot -100 to 400 for an epoch of 300. There are large shifts at the beginning and end (although ramp ups should be dealt with in the data prep)
        tbandForPlot=tBand(bandPlotIdx);
        tBandHold=tbandForPlot-epoch-holdPre; %300ms before hover period ended 
        tBandMove=tbandForPlot+movePre; %200ms after the target appears
        
        
        iti=specTempHoldI.dataSpec.data(:, 1:size(holdDataAllTemp.dataSpec.data,2), :); %match the iti side

    
    %Plot the channels    
    idxCC=1;
    
    figtitle=['Smoothed spike mean '];
    figure('Name', figtitle, 'Position', [5 150 1200 750]) %x bottom left, y bottom left, x width, y height
    set(gca,'FontSize', 22)
    sgtitle(figtitle)
    for jj=1:length(ChSpk)
        chstr=['ch', num2str(ChSpk(jj))];
        subplot(2,1,1)
        hold on
        H2=plot(tspkMove, perCh.(chstr).meansteSpkMove(:,1));
        H2.LineWidth=2;
        H2.Color=C(idxCC,:);
        H2.DisplayName=chstr;
        title('Mean Smoothed Spike Rates Move by channel')
        legend
        
        subplot(2,1,2)
        hold on
        H2=plot(tspkHold, perCh.(chstr).meansteSpkHold(:,1));
        H2.LineWidth=2;
        H2.Color=C(idxCC,:);
        H2.DisplayName=chstr;
        title('Mean Smoothed Spike Rates Hold by channel')
        legend
        
        idxCC=idxCC+4;
        
    end
    

    %% bands
    for jj=1:length(ChSpk)
        chstr=['ch', num2str(ChSpk(jj))];%%
        figtitle=['Band filtered data ', chstr];
        figure('Name', figtitle, 'Position', [5 150 1200 750]) %x bottom left, y bottom left, x width, y height
        set(gca,'FontSize', 22)
        sgtitle(figtitle)
        
        %lblPhase={'Move'; 'Hold'; 'Transition to Move'; 'Transition to Hold'};
        lblPhase={'Move'; 'Hold'};

        for rr=1:4
            idxC=1;            
            subplot(2, 2, rr)            
            hold on
            H1=shadedErrorBar(tbandForPlot, perCh.(chstr).meansteBandMove(bandPlotIdx, rr), perCh.(chstr).meansteBandMove(bandPlotIdx, rr+4));
            H1.mainLine.LineWidth=2;
            H1.patch.FaceColor=C(colorShBand(idxC),:);
            H1.patch.EdgeColor=C(colorShBand(idxC)-1,:);
            H1.mainLine.Color=C(colorShBand(idxC)-1,:);
            H1.edge(1).Color=C(colorShBand(idxC)-1,:);
            H1.edge(2).Color=C(colorShBand(idxC)-1,:);
            H1.mainLine.DisplayName=lblPhase{1};
            title(lblB{rr})
            idxC=idxC+1;
            
            H2=shadedErrorBar(tbandForPlot, perCh.(chstr).meansteBandHold(bandPlotIdx, rr), perCh.(chstr).meansteBandHold(bandPlotIdx, rr+4));
            H2.mainLine.LineWidth=2;
            H2.patch.FaceColor=C(colorShBand(idxC),:);
            H2.patch.EdgeColor=C(colorShBand(idxC)-1,:);
            H2.mainLine.Color=C(colorShBand(idxC),:);
            H2.edge(1).Color=C(colorShBand(idxC)-1,:);
            H2.edge(2).Color=C(colorShBand(idxC)-1,:);
            H2.mainLine.DisplayName=lblPhase{2};
            ax=gca;
            ax.XLim=[floor(tbandForPlot(1)) ceil(tbandForPlot(end))];
            
            
            idxC=idxC+1;            
            if rr==2
                legend([H1.mainLine H2.mainLine ],{'Move', 'Hold'});
            end
            if rr==5
                xlabel('Time (ms)');
                ylabel('Arbitrary units');
            end
        end
    end
        
        
   
    for jj=1:length(ChSpk)
        idxC=1;
        chstr=['ch', num2str(ChSpk(jj))];
        figtitle=['Smoothed Spike Mean and SE ', chstr, ];
        figure('Name', figtitle, 'Position', [5 150 1200 750]) %x bottom left, y bottom left, x width, y height
        set(gca,'FontSize', 22)
        title(figtitle)
        hold on
        H1=shadedErrorBar(tSpk, perCh.(chstr).meansteSpkMove(:,1), perCh.(chstr).meansteSpkMove(:,2));
        H1.mainLine.LineWidth=4;
        H1.patch.FaceColor=C(colorSh(idxC),:);
        H1.patch.EdgeColor=C(colorSh(idxC)-1,:);
        H1.mainLine.Color=C(colorSh(idxC)-1,:);
        H1.edge(1).Color=C(colorSh(idxC)-1,:);
        H1.edge(2).Color=C(colorSh(idxC)-1,:);
        H1.mainLine.DisplayName='Move';
        idxC=idxC+1;
        
        H2=shadedErrorBar(tSpk, perCh.(chstr).meansteSpkHold(:,1), perCh.(chstr).meansteSpkHold(:,2));
        H2.mainLine.LineWidth=4;
        H2.patch.FaceColor=C(colorSh(idxC),:);
        H2.patch.EdgeColor=C(colorSh(idxC)-1,:);
        H2.mainLine.Color=C(colorSh(idxC)-1,:);
        H2.edge(1).Color=C(colorSh(idxC)-1,:);
        H2.edge(2).Color=C(colorSh(idxC)-1,:);
        H2.mainLine.DisplayName='Hold';

        idxC=idxC+1;
        legend([H1.mainLine H2.mainLine ],{'Move', 'Hold'});

    end  
    %%
    for jj=1:length(ChSpk)
        
        chstr=['ch', num2str(ChSpk(jj))];
        %set up easy plotting
        Hh=perCh.(chstr).specHold;
        Mm=perCh.(chstr).specMove;
        perCh.(chstr).HoldMove=Hh-Mm;
        perCh.(chstr).HoldIti=Hh-iti(:,:,ChSpk(jj));
        perCh.(chstr).MoveIti=Mm-iti(:,:,ChSpk(jj));
        MmHh=Mm-Hh;
        HhI=Hh-iti(:,:,ChSpk(jj));
        MmI=Mm-iti(:,:,ChSpk(jj));
        %NEED TO MAKE THESE LOOK NICER, LIKELY REMOVE THE iti ONE OR THE NO
        %ITI ONE
       
        chstr=['ch', num2str(ChSpk(jj))];
        figtitle=['Difference between epochs ', chstr, ' Heatmap ', ];
        figure('Name', figtitle, 'Position', [5 150 1200 750]) %x bottom left, y bottom left, x width, y height
        set(gca,'FontSize', 22)
        sgtitle(figtitle)
        
%         subplot(3,2,2)
%         imagesc(tpHold, fDiff, HhI(:, specPlotIdx)); axis xy;
%         title('Hold - ITI')
%         colorbar
%         ax=gca;
%         ax.YTick=(0:40:round(fDiff(end)));
        
%         subplot(3,2,4)
%         imagesc(tpMove, fDiff, MmI(:, specPlotIdx)); axis xy;
%         title('Move - ITI')
%         colorbar
%         ax=gca;
%         ax.YTick=(0:40:round(fDiff(end)));
%         

        subplot(3, 1 ,1)
        imagesc(tpHold, fDiff, Hh(:, specPlotIdx)); axis xy;
        title('Hold')
        colorbar
        ax=gca;
        ax.YTick=(0:40:round(fDiff(end)));
        
        subplot(2, 1, 2)
        imagesc(tpMove, fDiff, Mm(:, specPlotIdx)); axis xy;
        title('Move')
        colorbar
        ax=gca;
        ax.YTick=(0:40:round(fDiff(end)));
        
        subplot(3,1, 3)
        ax=gca;
        imagesc(tpDiffforPlot, fDiff, MmHh(:, specPlotIdx)); axis xy;
        title('Move-Hold')
        colorbar
        ax=gca;
        ax.YTick=(0:40:round(fDiff(end)));

    end
end

%%

    


end

