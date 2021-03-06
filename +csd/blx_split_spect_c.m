function [ Ss,ts, fs ] = blx_split_spect_c( data, params, varargin )
%blx_split_spect_c Sets up a moving mtspectgram to take at regular time
%intervals and to subtract a baseline.
% %   INPUTS
%         data=   the data type, can be anything
%         map=    the map for titles
%     Varargins
%         ch=     the channels you want to run through this
%         t=      the time you want to run through
%     OUTPUTS
%         Ss=     an array of S outputs from spectrumc
%         fs=     an array of f outputs from spectrumc
% 
%
% %EXAMPLE:
% params=[];
% ch1n=map.ChannelInfo{ch1,2}; ch2n=map.ChannelInfo{ch2,2};
% figtitle=[blc.SourceBasename ' Spectrogram as compared to first 5 min baseline ' ch1n{1,1} ' to ' ch2n{1,1}]; %alternatively 'Bridget Miller_Sz2_6212016_Ph2D6_005855-000 Spectrogram as a ratio compared to first 5 min baseline PST1 to PST6'
%  [Ss,ts, fs ] = csd.blx_split_spect_c( data, params, 'tt', tt11PST, 'gridtype', 2, 'figtitle', figtitle, 'map', map )
% %
%   WORKS AS OF 8/29/17.  MAY NEED TO BE TESTED WITH A SET THAT HAS MORE
%   VARIABILITY (THE TEST WAS PRETTY UNIFORM SO HAD TO TELL IF IT WAS THE
%   ANALYSIS OF THE SUBJECT
%%
if nargin == 1 || isempty(params) % no parameter structure provided
    params = struct;
    params.Fs = 400;   % in Sampling rate (not hz)
    params.fpass = [0 70];     % [minFreqz maxFreq] in Hz
    params.tapers = [5 9];
    params.pad = 0;
    params.trialave = 0;
    params.win = [5 5];% size and step size for windowing continuous data, in seconds.  if you make the step size the same as the window size it won't skip any data.
    params.baseline=300; %size of the baseline window you want, in seconds (so 300 s is 5 min)
end



%%
%check if number of channels was specified
[varargin, ch]=util.argkeyval('ch', varargin, [1 size(data,2)]);
ch1=ch(1,1);
ch2=ch(1,2);
%%
%check if time plot was specified
[varargin, tt]=util.argkeyval('tt', varargin, []);
%%
%check if number of channels was specified, 1 is a 4x5 and 2 is a 1x6
[varargin, gridtype]=util.argkeyval('gridtype', varargin, 1);
%%
%check if the title was specified
[varargin, figtitle, ~, found ]=util.argkeyval('figtitle', varargin, []);
if ~found
    figtitle=['S0? Sz? spectrogram as ratio of baseline at first 5 min'];
end
%%
%make the channel names from the map
%check if the map was specified
[varargin, map, ~, found ]=util.argkeyval('map', varargin, []);
if ~found
    chname=ch(1,1):ch(1,2)+1;
    chname=arrayfun(@num2str,chname,'uniformoutput',false);
else
    chname=map.ChannelInfo{ch(1,1):ch(1,2)+1,2};
end
%%
switch gridtype
    case 1 %4x5 grid
        figure('Name', figtitle,'NumberTitle', 'off')
        idx=1;
        for ii=1:size(data,2)
            if ii==5 || ii== 10 || ii==15
                continue
            else
                subplot(4,4,idx)
                if isempty(tt)
                    [S,tt,fs] = mtspecgramc(data(:,ii), params.win, params);
                    %can change this to whatever, remember that the window
                    %is in seconds. This does a single spectrogram chunk
                    [Ss]=mtspecgramc(data(1:params.baseline*400, ii), [params.baseline params.baseline], params);
                    for jj=1:size(S,2)
                        % makes S a ratio of the baseline
                        S(:,jj)=S(:,jj)/Ss(1,jj);
                    end
                else
                    [S,ts,fs] = mtspecgramc(data(:,ii), params.win, params);
                    [Ss]=mtspecgramc(data(1:params.baseline*400, ii), [params.baseline params.baseline], params);
                    for jj=1:size(S,2)
                        % makes S a ratio of the baseline
                        S(:,jj)=S(:,jj)/Ss(1,jj);
                    end
                end
                imagesc(tt,fs, 10*log10(S)');
                axis xy;
                
                title([' Ch ' chname{ii,1} '-' chname{ii+1,1}]);
                idx=idx+1;
            end
        end
    case 2 %6 electrode strip
        figure('Name', figtitle,'NumberTitle', 'off')
        idx=1;
        for ii=1:size(data,2)
            subplot(5,1,idx)
            if isempty(tt)
                [S,tt,fs] = mtspecgramc(data(:,ii), params.win, params);
                [Ss]=mtspecgramc(data(1:params.baseline*400, ii), [params.baseline params.baseline], params);
                for jj=1:size(S,2)
                    % makes S a ratio of the baseline
                    S(:,jj)=S(:,jj)/Ss(1,jj);
                end
            else
                [S,ts,fs] = mtspecgramc(data.data(:,ii), params.win, params);
                [Ss]=mtspecgramc(data(1:params.baseline*400, ii), [params.baseline params.baseline], params);
                for jj=1:size(S,2)
                    % makes S a ratio of the baseline
                    S(:,jj)=S(:,jj)/Ss(1,jj);
                end
            end
            imagesc(tt,fs, 10*log10(S)');
            axis xy;
            
            title([' Ch ' chname{ii,1} '-' chname{ii+1,1}]);
            idx=idx+1;
        end
        case 3 %3 1x6 electrode strips
        figure('Name', figtitle,'NumberTitle', 'off')
        idx=1;
        for ii=1:size(data,2)
            if ii==6 || ii== 12              
                continue
            else
                subplot(3,5,idx)
                if isempty(tt)
                    [S,tt,fs] = mtspecgramc(data(:,ii), params.win, params);
                    [Ss]=mtspecgramc(data(1:params.baseline*400, ii), [params.baseline params.baseline], params);
                    for jj=1:size(S,2)
                        % makes S a ratio of the baseline
                        S(:,jj)=S(:,jj)/Ss(1,jj);
                    end
                else
                    [S,ts,fs] = mtspecgramc(data(:,ii), params.win, params);
                    [Ss]=mtspecgramc(data(1:params.baseline*400, ii), [params.baseline params.baseline], params);
                    for jj=1:size(S,2)
                        % makes S a ratio of the baseline
                        S(:,jj)=S(:,jj)/Ss(1,jj);
                    end
                end
                imagesc(tt,fs, 10*log10(S)');
                axis xy;
                
                title([' Ch ' chname{ii,1} '-' chname{ii+1,1}]);
                idx=idx+1;
            end
        end
        
end



end



