function plot_spec_target_avg(specs, timebins, freqbins, channels, channel_key, varargin)
    
    [varargin, title_mod, ~, ~] = util.argkeyval('TitleMod', varargin, '');
    [varargin, text_add, ~, add_text] = util.argkeyval('AddText', varargin, '');
    [varargin, freq_range, ~, freq_range_found] = util.argkeyval('FreqRange', varargin, []);
    util.argempty(varargin);
    
    if freq_range_found
        freq_indx = freqbins >= freq_range(1) & freqbins <= freq_range(2);
        freqbins = freqbins(freq_indx);
        specs = specs(:, freq_indx, :, :);
        title_mod = sprintf('%s %d-%dHz',title_mod, freq_range(1), freq_range(2));
        
    end
    
    num_chans = size(specs, 3);
    assert(num_chans == length(channels), 'channels dont match my friend')
    t_subplot_order = [8 1 2 7 0 3 6 5 4];
    
    
    
    for ch = 1:num_chans
        figtitle = sprintf('All Targs%s - ch-%d - %s', title_mod, channels(ch), channel_key(ch));
        figure('Name', figtitle, 'NumberTitle', 'off','position', [-1919 121 1920 1083])
        
        for t = 1:length(t_subplot_order)
            if t == 5
                continue %skip the middle subplot
            end
            target = t_subplot_order(t);
            subplot_tight(3, 3, t, [0.08 0.08])
            imagesc(timebins, freqbins, specs(:, :, ch, target)'); axis xy;
            ylabel('Frequency (Hz)')
            xlabel('Time (s)')
            ts = sprintf('Ch %d S %d Target %d', channels(ch), ch, target);
            title(ts)
            colorbar
        end % end for subplots
        if add_text
            annotation('textbox', [0.40544 0.46831 0.18705 0.043062], 'String', text_add, 'FitBoxToText', 'on')
        end
    end % end for channels
end % end func get_avg_and_diff