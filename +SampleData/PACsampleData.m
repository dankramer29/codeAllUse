function [signalfinal] = PACsampleData(freqRange, lowFreq, highFreq)

%Creates sample Phase Amplitude Coupling data.
% 
% Inputs:
%     freqRange= [startingfreq endfreq] typically 1 200 will work
%     lowFreq= the low frequency you want to couple, can be a band
%     highFreq= the high frequency you want to couple, can be a band
% Outputs 
%    signalfinal= the output of a single signal with noise built in,
% and coupling of the two frequencies
% 
% couples by multiplying the higher frequency by increasing factors during
% only the up phase of the lower frequency

%to create sample data form this for ucsfPAC, run this:
%{
for kk=1:8
    for jj=1:9
        signal(:,jj,kk)=Analysis.SampleData.PACsampleData(1:200, 20, 70);
    end
end
%}


%TO DO, MAKE IT SO YOU AN SEE A BAND OF FREQUENCIES

combined=true;

plotOn=false; %turn to true if you want to plot the fake data

fs=1/2000; %sampling rate
T=3; %lenth of recording time, in seconds
t=0:fs:T; %x input, i.e. the points

%how much to divide the amplitude of the noise by, so that it's not a crazy
%amount, if 2, you are dividing the amplitude of the noise by half the
%amplitude of that frequency
noiseFactor=2; 


%% SET FREQUENCY PARAMETERS
%frequency of your signal
%complex signal,builds a signal from your low to high frequencies, starting
%at 1Hz
f=freqRange; %can also make independent fs, e.g. f=[20 80]
AmpBase=0.25; %the amplitude range you want
signalfinal=zeros(1,length(t));
idx1=1;
for ii=1:length(f)
    clear ampF, clear xx, clear zz, clear signalfreq, clear noise,
    f(ii)=ii;
%     if ii==1 || ii==2 || ii==3
%         Amp(ii)=AmpBase*.1;
%     else
        Amp(ii)=AmpBase*1/f(ii); %amplitude, 1/f
%     end
    noise=Amp(ii)/noiseFactor*(rand(1,length(t)));
    phse=randi(360);
    
    if ismember(ii,lowFreq)
        phse=0;
    end
    if ismember(ii,highFreq)
         
        phse=0;
        halfphase=0:fs:1/(lowFreq(idx1)*2)-fs; %create time points for half the phase of the lower freq
        xx=sin(2*pi*lowFreq(idx1)*halfphase); %make a sin wave for half the phase
        zz=zeros(1,length(halfphase)); %make ones for the other half
        xx=[xx zz]; %combine them 
        ampF=repmat(xx,1,lowFreq(idx1)*T); %make this the length of the time series
        dif=length(t)-length(ampF); %whatever is short from the length, add zeros, (e.g. 21 Hz isn't divisible into 2000 samples/sec)
        ampF=[ampF zeros(1,dif)];
        ampF=ampF+1; %add 1 so that this becomes a multiplication factor to amplify only the in phase stuff
       
        signalfreq=Amp(ii)*sin(2*pi*f(ii)*t+phse);
        signalfreq=signalfreq+noise; %add the noise
        signalfreq=signalfreq.*ampF;
        signalPh(:,idx1)=signalfreq;
        idx1=idx1+1;
        
    else
        signalfreq=Amp(ii)*sin(2*pi*f(ii)*t+phse);
        signalfreq=signalfreq+noise; %add the noise
    end

    
    signalfinal=signalfinal+signalfreq;
end



%%
if plotOn
    figure
    hold on
    plot(t,signalfinal)
end
