%% THIS IS TO OPEN THE FILE, AND RUN IT AS A DATAPROC SO YOU DON'T HAVE TO KEEP RUNNING IT EACH TIME FOR DEBUGGING
%TO RUN THESE INDIVIDUALLY, HIT RUN SECTION AND YOU CAN DEBUG

%Task open YES USE THIS ONE, IT'S THE SERGEY ONE
ns5file='C:\Users\dankr\Documents\Data\t5.2019.03.18_HeadMoveData\Data\NS5_block10Head\10_cursorTask_Complete_t5_bld(010).ns5';
load('C:\Users\dankr\Documents\Data\t5.2019.03.18_HeadMoveData\Data\R_structs\R_block_10.mat')
R=R(1:64);
%optional
ns5=openNSx('read', ns5file, 'c:01:96'); 
data=double(ns5.Data{end})'; %may need {3} depending on the data type


%%
R=R(1:10);
rState=double([R.state]);
rCP=double([R.cursorPosition]); %1 and 3 are position, 2 and 4 are velocity for grid task, and it's derivative of cp in grid task or row 1:2 of xp
rFCT=double([R.firstCerebusTime])+ 10*30; %There is a 10ms offset from ns5 and firstCerebusTime, unknown why, but verified. difference between row 1 and 2 is 6ms, lastcerebrustime is 29 and 35ms respectively for row 1/2
rMinA=double([R.minAcausSpikeBand])';
if isempty(cell2mat(strfind(fields(R), 'spikeRaster')))
    threshMultiplier = -4.5;
    RMSvals = channelRMS( R );
    thresholds = RMSvals.*threshMultiplier;
    thresholds = thresholds'; % make column vector
    for iTrial = 1 : numel( R )
        R(iTrial).spikeRaster = R(iTrial).minAcausSpikeBand < repmat( thresholds, 1, size( R(iTrial).minAcausSpikeBand, 2 ) );
    end
end
rSpk=double([R.spikeRaster]);

%% prep all the data.
fs=30000;
preTrialSt=2; %start 2 seconds before the trial
preTrialEnd=2; %start 2 seconds after the trial
trialSt=rFCT(1)-fs*preTrialSt;
trialEnd=rFCT(end)+fs*preTrialEnd;
downSample=2000;
tT=tic;
[dataProc.DataAllTemp, dataProc.params]=Analysis.BasicDataProc.dataPrep(Analysis.BasicDataProc.downSample(data, fs, downSample), 'Spectrogram', true, 'doBandFilterBroad', true);
toc(tT)