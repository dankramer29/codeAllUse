function plotImpedanceCentral(impedanceFile,mapFile)
% plots impedances like Central does

hArray = Blackrock.ArrayMap(mapFile);
[impedance,channels,electrodes] = Blackrock.processImpedance(impedanceFile,mapFile);

h = figure('Name',impedanceFile,...
    'PaperPositionMode','auto',...
    'Position',[50 50 1200 800]);
sph = hArray.getChannelSubplot(channels,'Parent',h);

for k=channels(:)'
    axis(sph(k),'on');
    if impedance(k) > 800
        set(sph(k),'Color',[1 0 0]);
        textColor = [0.9 0.9 0.9];
    elseif impedance(k) <= 50
        set(sph(k),'Color',[1 1 0.5]);
        textColor = [0.1 0.1 0.1];
    else
        set(sph(k),'Color',[0 0 0]);
        textColor = [0.9 0.9 0.9];
    end
    set(sph(k),'XTick',[],'YTick',[],'Box','on');
    text(0.5,0.5,{num2str(impedance(k)),'kOhm'},...
        'Parent',sph(k),...
        'FontSize',10,...
        'FontWeight','bold',...
        'HorizontalAlignment','center',...
        'Color',textColor);
    text(0.02,0.92,['chan' num2str(channels(k))],...
        'Parent',sph(k),...
        'FontSize',8,...
        'FontWeight','normal',...
        'Color',textColor);
    text(0.02,0.78,['elec' num2str(electrodes(k))],...
        'Parent',sph(k),...
        'FontSize',8,...
        'FontWeight','normal',...
        'Color',textColor);
end
