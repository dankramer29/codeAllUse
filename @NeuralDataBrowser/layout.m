function layout(this)

% delete any old GUIs
ff = findobj('Name',this.name);
delete(ff);

% set figure position based on screen dimensions
figpos = [250 40 this.width this.height];
rootProps = get(0);
if isfield(rootProps,'ScreenSize')
    figleft = max(rootProps.ScreenSize(1)+this.screenMargin(1)-1,(rootProps.ScreenSize(3)-this.width-this.screenMargin(3))/2);
    figbottom = max(rootProps.ScreenSize(2)+this.screenMargin(2)-1,(rootProps.ScreenSize(4)-this.height-this.screenMargin(4))/2);
    figpos = [figleft figbottom this.width this.height];
end

this.axisHeight = round(0.6*(this.height - this.titleHeight - 2*this.axisSpacing - 2*this.outerSpacing));
configHeight = this.height - this.axisHeight - this.axisConfigSpacing - this.titleHeight - 2*this.axisSpacing - 2*this.outerSpacing;

% create the figure
this.hFigure = figure(...
    'Units','pixels',...
    'Color',[0.94 0.94 0.94],...
    'Position',figpos,...
    'PaperPositionMode','auto',...
    'NumberTitle','off',...
    'WindowKeyPressFcn',@(h,kp)keyPressHandler(this,kp),...
    'Resize','off',...
    'MenuBar','none',...
    'name',this.name,...
    'ToolBar','none',...
    'Tag','ndb');
this.hFigure.CloseRequestFcn = @(src,dt)close(this);

% status readout
currLeft = 10;
currBottom = 0;
localWidth = this.width-2*10;
localHeight = 25;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom localWidth localHeight],...
    'Enable','off',...
    'String','',...
    'HorizontalAlignment','left',...
    'Style','text',...
    'FontSize',10,...
    'Tag','textStatus');

% neural data directory
currLeft = this.outerSpacing;
currBottom = this.outerSpacing + configHeight + this.axisConfigSpacing + this.axisSpacing + this.axisHeight + this.axisSpacing;
localWidth = 50;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','left',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Data Dir: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 200;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','left',...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String','[No Directory Selected]',...
    'Style','edit',...
    'Tag','editDataDirectory');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 24;
localHeight = localWidth;
classdir = fileparts(mfilename('fullpath'));
imgdata = imread(fullfile(classdir,'img/folder-open-icon-24.png'));
grayval = get(gcf,'Color');
for kk=1:size(imgdata,3)
    tmp = imgdata(:,:,kk);
    tmp(tmp==0)=round(grayval(kk)*255);
    imgdata(:,:,kk)=tmp;
end
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonSelectDataFolder',...
    'CData',imgdata,...
    'Callback',@(h,evt)buttonSelectDataFolder_Callback);

% neural data files
currLeft = currLeft + localWidth + 6*this.elemSpacing;
localWidth = 55;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Data File: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = this.width - this.outerSpacing - this.elemSpacing - currLeft - (this.elemSpacing+24);
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',{'No Files Available'},...
    'Enable','off',...
    'Style','popupmenu',...
    'Tag','popupDataFiles',...
    'Callback',@(h,evt)true);
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 24;
localHeight = localWidth;
classdir = fileparts(mfilename('fullpath'));
imgdata = imread(fullfile(classdir,'img/Actions-dialog-ok-apply-icon-24.png'));
grayval = get(gcf,'Color');
for kk=1:size(imgdata,3)
    tmp = imgdata(:,:,kk);
    tmp(tmp==0)=round(grayval(kk)*255);
    imgdata(:,:,kk)=tmp;
end
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonLoadFile',...
    'CData',imgdata,...
    'Callback',@(h,evt)buttonLoadDataFile_Callback);








% task directory
currLeft = this.outerSpacing;
currBottom = currBottom - this.rowHeight - 2;
localWidth = 50;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','left',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Task Dir: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 200;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','left',...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String','[No Directory Selected]',...
    'Style','edit',...
    'Tag','editTaskDirectory');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 24;
localHeight = localWidth;
classdir = fileparts(mfilename('fullpath'));
imgdata = imread(fullfile(classdir,'img/folder-open-icon-24.png'));
grayval = get(gcf,'Color');
for kk=1:size(imgdata,3)
    tmp = imgdata(:,:,kk);
    tmp(tmp==0)=round(grayval(kk)*255);
    imgdata(:,:,kk)=tmp;
end
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonSelectTaskFolder',...
    'CData',imgdata,...
    'Callback',@(h,evt)buttonSelectTaskFolder_Callback);

% task files
currLeft = currLeft + localWidth + 6*this.elemSpacing;
localWidth = 55;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Task File: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = this.width - this.outerSpacing - this.elemSpacing - currLeft - (this.elemSpacing+24);
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',{'No Files Available'},...
    'Enable','off',...
    'Style','popupmenu',...
    'Tag','popupTaskFiles',...
    'Callback',@(h,evt)true);
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 24;
localHeight = localWidth;
classdir = fileparts(mfilename('fullpath'));
imgdata = imread(fullfile(classdir,'img/Actions-dialog-ok-apply-icon-24.png'));
grayval = get(gcf,'Color');
for kk=1:size(imgdata,3)
    tmp = imgdata(:,:,kk);
    tmp(tmp==0)=round(grayval(kk)*255);
    imgdata(:,:,kk)=tmp;
end
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonLoadTaskFile',...
    'CData',imgdata,...
    'Callback',@(h,evt)buttonLoadTaskFile_Callback);

popupTaskPos = get(findall(this.hFigure,'tag','popupTaskFiles'),'position');
currLeft = popupTaskPos(1);
localHeight = 14;
currBottom = popupTaskPos(2)-localHeight;
localWidth = popupTaskPos(3);
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','left',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','',...
    'Tag','textDataInfo',...
    'Style','text');




% main axis
textInfoPos = get(findall(this.hFigure,'tag','textDataInfo'),'position');
currLeft = this.elemSpacing + this.axisLeftSpacing;
currBottom = currBottom - this.axisHeight - this.axisSpacing + textInfoPos(4);
localWidth = this.width - 2*this.elemSpacing - this.axisLeftSpacing - this.outerSpacing;
localHeight = this.axisHeight;
axes(...
    'Parent',this.hFigure,...
    'Units','pixels',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'Box','on',...
    'Tag','axisMain');

% navigation
navBtnWidth = 100;
navBtnHeight = 3*this.rowHeight;

% backward
currLeft = this.elemSpacing + this.axisLeftSpacing;
axisPos = get(findall(this.hFigure,'tag','axisMain'),'Position');
currBottom = axisPos(2) - this.axisConfigSpacing - navBtnHeight;
localWidth = navBtnWidth;
localHeight = navBtnHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonNavigateBeginning',...
    'String',' |< ',...
    'FontSize',16,...
    'Callback',@(h,evt)buttonNavigate_Callback(-inf));
currLeft = currLeft + localWidth + this.elemSpacing;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonNavigateBackward',...
    'String',' < ',...
    'FontSize',16,...
    'Callback',@(h,evt)buttonNavigate_Callback(-this.timeStep));

% forward
currLeft = this.width - this.outerSpacing - 2*localWidth - this.elemSpacing;
localWidth = navBtnWidth;
localHeight = navBtnHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonNavigateForward',...
    'String',' > ',...
    'FontSize',16,...
    'Callback',@(h,evt)buttonNavigate_Callback(this.timeStep));
currLeft = currLeft + localWidth + this.elemSpacing;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+2 localWidth localHeight],...
    'Style','pushbutton',...
    'Tag','buttonNavigateEnd',...
    'String',' >| ',...
    'FontSize',16,...
    'Callback',@(h,evt)buttonNavigate_Callback(inf));

% time range
navBegPos = get(findall(this.hFigure,'tag','buttonNavigateBeginning'),'Position');
currLeft = this.elemSpacing + this.axisLeftSpacing;
currBottom = navBegPos(2) - this.rowHeight - 5*this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Time Range: ',...
    'Style','text');
currLeft = currLeft + navBtnWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.timeRange),...
    'Style','edit',...
    'Tag','editTimeRange',...
    'Callback',@(h,evt)editTimeRange_Callback,...
    'KeypressFcn',@(h,evt)editTimeRange_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyTimeRange',...
    'Callback',@(h,evt)buttonApplyTimeRange_Callback);

% time step
currLeft = this.elemSpacing + this.axisLeftSpacing;
currBottom = currBottom - this.rowHeight - this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Time Step: ',...
    'Style','text');
currLeft = currLeft + navBtnWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.timeStep),...
    'Style','edit',...
    'Tag','editTimeStep',...
    'Callback',@(h,evt)editTimeStep_Callback,...
    'KeypressFcn',@(h,evt)editTimeStep_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyTimeStep',...
    'Callback',@(h,evt)buttonApplyTimeStep_Callback);

% y limit min
currLeft = this.elemSpacing + this.axisLeftSpacing;
currBottom = currBottom - this.rowHeight - 3*this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','y-axis Minimum: ',...
    'Style','text');
currLeft = currLeft + navBtnWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.yLimitMinimum),...
    'Style','edit',...
    'Tag','editYLimitMinimum',...
    'Callback',@(h,evt)editYLimitMinimum_Callback,...
    'KeypressFcn',@(h,evt)editYLimitMinimum_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyYLimitMinimum',...
    'Callback',@(h,evt)buttonApplyYLimitMinimum_Callback);

% y limit max
currLeft = this.elemSpacing + this.axisLeftSpacing;
currBottom = currBottom - this.rowHeight - this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','y-axis Maximum: ',...
    'Style','text');
currLeft = currLeft + navBtnWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.yLimitMaximum),...
    'Style','edit',...
    'Tag','editYLimitMaximum',...
    'Callback',@(h,evt)editYLimitMaximum_Callback,...
    'KeypressFcn',@(h,evt)editYLimitMaximum_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyYLimitMaximum',...
    'Callback',@(h,evt)buttonApplyYLimitMaximum_Callback);

% y-limit auto/manual
currLeft = this.elemSpacing + this.axisLeftSpacing + 100 + this.elemSpacing;
currBottom = currBottom - this.rowHeight;
localHeight = this.rowHeight;
localWidth = 100;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[0.94 0.94 0.94],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String','Auto-Scale',...
    'Style','checkbox',...
    'Tag','checkboxYLimitAutoScale',...
    'Callback',@(h,evt)checkboxYLimitAutoScale_Callback);



% select all / none
chLayoutWidth = 350;
chLayoutHeight = 200;
navBackPos = get(findall(this.hFigure,'tag','buttonNavigateBackward'),'position');
navForwardPos = get(findall(this.hFigure,'tag','buttonNavigateForward'),'position');
statPos = get(findall(this.hFigure,'tag','textStatus'),'position');
betweenLR = navForwardPos(1) - (navBackPos(1)+navBackPos(3));
marginLR = (betweenLR - chLayoutWidth)/2;
currLeft = navBackPos(1)+navBackPos(3)+marginLR;
currBottom = statPos(2)+statPos(4) + 2;
localWidth = 174.5;
localHeight = 25;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Select All Channels',...
    'Style','pushbutton',...
    'Tag','buttonSelectAllChannels',...
    'Callback',@(h,evt)buttonSelectAllChannels_Callback);
currLeft = currLeft + localWidth + chLayoutWidth-2*localWidth;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Deselect All Channels',...
    'Style','pushbutton',...
    'Tag','buttonDeselectAllChannels',...
    'Callback',@(h,evt)buttonDeselectAllChannels_Callback);

% channel layout
selAllPos = get(findall(this.hFigure,'tag','buttonSelectAllChannels'),'position');
betweenUD = (navForwardPos(2)+navForwardPos(4)-50) - (selAllPos(2)+selAllPos(4)+2);
marginUD = (betweenUD - chLayoutHeight)/2;
currLeft = navBackPos(1)+navBackPos(3)+marginLR;
currBottom = selAllPos(2)+selAllPos(4)+marginUD;
chSpacing = 1;
chWidth = chLayoutWidth/10 - chSpacing;
chHeight = chLayoutHeight/10 - chSpacing;
layout = this.hArrayMap.vec2layout(1:96);
for kk=1:96
    [row,col] = find(layout==kk);
    btnpos = [currLeft + (col-1)*(chWidth+chSpacing) currBottom + (10-row+1)*(chHeight+chSpacing) chWidth chHeight];
    uicontrol(...
        'Parent',this.hFigure,...
        'Position',btnpos,...
        'Enable','off',...
        'String',kk,...
        'Style','togglebutton',...
        'Tag',sprintf('button%02d',kk),...
        'Callback',@(h,evt)buttonChannel_Callback(kk));
end

% plot type
navForwardPos = get(findall(this.hFigure,'tag','buttonNavigateForward'),'position');
currLeft = navForwardPos(1);
currBottom = navForwardPos(2) - this.rowHeight - 5*this.elemSpacing;
localWidth = 48;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Type: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 100;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',{'broadband','mn+std','distributed','spectrogram'},...
    'Style','popupmenu',...
    'Tag','popupPlotType',...
    'Callback',@(h,evt)popupPlotType_Callback);
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyPlotType',...
    'Callback',@(h,evt)buttonApplyPlotType_Callback);

% window
currLeft = this.width - this.outerSpacing - 2*navBtnWidth - this.elemSpacing;
currBottom = currBottom - this.rowHeight - 3*this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Window Size: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.window),...
    'Style','edit',...
    'Tag','editWindow',...
    'Callback',@(h,evt)editWindow_Callback,...
    'KeypressFcn',@(h,evt)editWindow_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyWindow',...
    'Callback',@(h,evt)buttonApplyWindow_Callback);

% step
currLeft = this.width - this.outerSpacing - 2*navBtnWidth - this.elemSpacing;
currBottom = currBottom - this.rowHeight - this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Step Size: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.step),...
    'Style','edit',...
    'Tag','editStep',...
    'Callback',@(h,evt)editStep_Callback,...
    'KeypressFcn',@(h,evt)editStep_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyStep',...
    'Callback',@(h,evt)buttonApplyStep_Callback);

% tapers
currLeft = this.width - this.outerSpacing - 2*navBtnWidth - this.elemSpacing;
currBottom = currBottom - this.rowHeight - this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Tapers: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',util.vec2str(this.tapers),...
    'Style','edit',...
    'Tag','editTapers',...
    'Callback',@(h,evt)editTapers_Callback,...
    'KeypressFcn',@(h,evt)editTapers_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyTapers',...
    'Callback',@(h,evt)buttonApplyTapers_Callback);

% pad
currLeft = this.width - this.outerSpacing - 2*navBtnWidth - this.elemSpacing;
currBottom = currBottom - this.rowHeight - this.elemSpacing;
localWidth = 100;
localHeight = this.rowHeight;
uicontrol(...
    'Parent',this.hFigure,...
    'HorizontalAlignment','right',...
    'Position',[currLeft currBottom localWidth localHeight],...
    'String','Zero-Padding: ',...
    'Style','text');
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 48;
uicontrol(...
    'Parent',this.hFigure,...
    'enable','off',...
    'HorizontalAlignment','right',...
    'BackgroundColor',[1 1 1],...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'String',num2str(this.pad),...
    'Style','edit',...
    'Tag','editPad',...
    'Callback',@(h,evt)editPad_Callback,...
    'KeypressFcn',@(h,evt)editPad_KeypressFcn(evt));
currLeft = currLeft + localWidth + this.elemSpacing;
localWidth = 47;
uicontrol(...
    'Parent',this.hFigure,...
    'Position',[currLeft currBottom+4 localWidth localHeight],...
    'Enable','off',...
    'String','Apply',...
    'Style','pushbutton',...
    'Tag','buttonApplyPad',...
    'Callback',@(h,evt)buttonApplyPad_Callback);


% pull out gui handles
this.guiHandles = guihandles(this.hFigure);



    function buttonSelectDataFolder_Callback
        datapath = env.get('data');
        try
            dirname = uigetdir(datapath{1},'Select Data Directory');
            if ischar(dirname)
                this.dataDirectory = dirname;
                updateGUI_ChangeDataDirectory(this);
            end
        catch ME
            util.errorMessage(ME);
        end
    end % END function buttonSelectFolder_Callback

    function buttonLoadDataFile_Callback
        loadDataFile(this);
    end % END function buttonLoadDataFile_Callback

    function buttonSelectTaskFolder_Callback
        datapath = env.get('data');
        try
            dirname = uigetdir(datapath{1},'Select Data Directory');
            if ischar(dirname)
                this.dataDirectory = dirname;
                updateGUI_ChangeTaskDirectory(this);
            end
        catch ME
            util.errorMessage(ME);
        end
    end % END function buttonSelectTaskFolder_Callback

    function buttonLoadTaskFile_Callback
        loadTaskFile(this);
    end % END function buttonLoadDataFile_Callback

    function buttonNavigate_Callback(amount)
        navigate(this,amount);
    end % END function buttonNavigateForward_Callback

    function editTimeRange_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTimeRange;
        hApply = this.guiHandles.buttonApplyTimeRange;
        eaEditCallback(hEdit,hApply,'timeRange',true);
    end % END function editTimeRange_Callback

    function editTimeRange_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTimeRange;
        eaEditKeypress(evt,hEdit);
    end % END function editTimeRange_KeypressFcn

    function buttonApplyTimeRange_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTimeRange;
        hApply = this.guiHandles.buttonApplyTimeRange;
        eaApplyCallback(hEdit,hApply,'timeRange',true);
    end % END function buttonApplyLimit_Callback

    function editTimeStep_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTimeStep;
        hApply = this.guiHandles.buttonApplyTimeStep;
        eaEditCallback(hEdit,hApply,'timeStep',true);
    end % END function editTimeStep_Callback

    function editTimeStep_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTimeStep;
        eaEditKeypress(evt,hEdit);
    end % END function editTimeStep_KeypressFcn

    function buttonApplyTimeStep_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTimeStep;
        hApply = this.guiHandles.buttonApplyTimeStep;
        eaApplyCallback(hEdit,hApply,'timeStep',true);
    end % END function buttonApplyTimeStep_Callback

    function editYLimitMinimum_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editYLimitMinimum;
        hApply = this.guiHandles.buttonApplyYLimitMinimum;
        eaEditCallback(hEdit,hApply,'yLimitMinimum',true);
    end % END function editYLimitMinimum_Callback

    function editYLimitMinimum_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editYLimitMinimum;
        eaEditKeypress(evt,hEdit);
    end % END function editYLimitMinimum_KeypressFcn

    function buttonApplyYLimitMinimum_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editYLimitMinimum;
        hApply = this.guiHandles.buttonApplyYLimitMinimum;
        eaApplyCallback(hEdit,hApply,'yLimitMinimum',true);
    end % END function buttonApplyYLimitMinimum_Callback

    function editYLimitMaximum_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editYLimitMaximum;
        hApply = this.guiHandles.buttonApplyYLimitMaximum;
        eaEditCallback(hEdit,hApply,'yLimitMaximum',true);
    end % END function editYLimitMaximum_Callback

    function editYLimitMaximum_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editYLimitMaximum;
        eaEditKeypress(evt,hEdit);
    end % END function editYLimitMaximum_KeypressFcn

    function buttonApplyYLimitMaximum_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editYLimitMaximum;
        hApply = this.guiHandles.buttonApplyYLimitMaximum;
        eaApplyCallback(hEdit,hApply,'yLimitMaximum',true);
    end % END function buttonApplyYLimitMaximum_Callback

    function checkboxYLimitAutoScale_Callback
        
        this.yLimitAutoScale = get(this.guiHandles.checkboxYLimitAutoScale,'value');
        updateGUI_DisplayFile(this);
    end % END checkboxYLimitAutoScale_Callback

    function buttonChannel_Callback(ch)
        
        oldval = get(this.guiHandles.(sprintf('button%02d',ch)),'Value');
        if isequal(oldval,get(this.guiHandles.(sprintf('button%02d',ch)),'Max'))
            
            % was not present, add
            this.currentChannels = unique([ch this.currentChannels]);
        else
            
            % was on, remove
            this.currentChannels = setdiff(this.currentChannels,ch);
        end
        updateGUI_DisplayFile(this);
    end % END function buttonChannel_Callback

    function buttonSelectAllChannels_Callback
        this.currentChannels = 1:96;
        updateGUI_DisplayFile(this);
    end % END function buttonSelectAllChannels_Callback

    function buttonDeselectAllChannels_Callback
        this.currentChannels = [];
        updateGUI_DisplayFile(this);
    end % END function buttonDeselectAllChannels_Callback

    function popupPlotType_Callback
        which = get(this.guiHandles.popupPlotType,'Value');
        list = get(this.guiHandles.popupPlotType,'String');
        val = list{which};
        if ~strcmpi(this.plotType,val)
            set(this.guiHandles.buttonApplyPlotType,'Enable','on');
        end
    end % END function popupPlotType_Callback

    function buttonApplyPlotType_Callback
        which = get(this.guiHandles.popupPlotType,'Value');
        list = get(this.guiHandles.popupPlotType,'String');
        val = list{which};
        this.plotType = val;
        set(this.guiHandles.buttonApplyPlotType,'Enable','off');
        updateGUI_DisplayFile(this);
    end % END buttonApplyPlotType_Callback

    function editWindow_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editWindow;
        hApply = this.guiHandles.buttonApplyWindow;
        eaEditCallback(hEdit,hApply,'window',true);
    end % END function editWindow_Callback

    function editWindow_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editWindow;
        eaEditKeypress(evt,hEdit);
    end % END function editWindow_KeypressFcn

    function buttonApplyWindow_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editWindow;
        hApply = this.guiHandles.buttonApplyWindow;
        eaApplyCallback(hEdit,hApply,'window',true);
    end % END function buttonApplyWindow_Callback

    function editStep_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editStep;
        hApply = this.guiHandles.buttonApplyStep;
        eaEditCallback(hEdit,hApply,'step',true);
    end % END function editStep_Callback

    function editStep_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editStep;
        eaEditKeypress(evt,hEdit);
    end % END function editStep_KeypressFcn

    function buttonApplyStep_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editStep;
        hApply = this.guiHandles.buttonApplyStep;
        eaApplyCallback(hEdit,hApply,'step',true);
    end % END function buttonApplyStep_Callback

    function editTapers_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTapers;
        hApply = this.guiHandles.buttonApplyTapers;
        eaEditCallback(hEdit,hApply,'tapers',true);
    end % END function editTapers_Callback

    function editTapers_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTapers;
        eaEditKeypress(evt,hEdit);
    end % END function editTapers_KeypressFcn

    function buttonApplyTapers_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editTapers;
        hApply = this.guiHandles.buttonApplyTapers;
        eaApplyCallback(hEdit,hApply,'tapers',true);
    end % END function buttonApplyTapers_Callback

    function editPad_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editPad;
        hApply = this.guiHandles.buttonApplyPad;
        eaEditCallback(hEdit,hApply,'pad',true);
    end % END function editPad_Callback

    function editPad_KeypressFcn(evt)
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editPad;
        eaEditKeypress(evt,hEdit);
    end % END function editPad_KeypressFcn

    function buttonApplyPad_Callback
        
        % pull together info and run the generic callback
        hEdit = this.guiHandles.editPad;
        hApply = this.guiHandles.buttonApplyPad;
        eaApplyCallback(hEdit,hApply,'pad',true);
    end % END function buttonApplyPad_Callback

    function eaEditCallback(hEdit,hApply,name,FlagNumeric)
        % EAEDITCALLBACK Generic edit callback for a specific context
        %
        %   EA stands for Edit / Apply.  It describes a specific setup of
        %   UI controls in which a user may change an edit box value, and
        %   either press enter or click an apply button.
        %
        %   This callback's purpose is to detect changes and either
        %   directly apply them (enter key pressed) or enable the apply
        %   button.
        
        % get new and old values
        newVal = get(hEdit,'String');
        if FlagNumeric
            if strcmpi(newVal(1),'[')
                newVal = eval(newVal);
            else
                newVal = str2double(newVal);
            end
        end
        oldVal = this.(name);
        
        
        % check if an updated value has been entered
        if (FlagNumeric && all(newVal == oldVal)) || (~FlagNumeric && strcmpi(newVal,oldVal))
            
            % same as old value, so no action
            set(hApply,'enable','off');
        else
            
            % detect enter/return, or enable apply button
            if get(hEdit,'Value')
                
                % keypress sets value=1, so directly update value
                set(hEdit,'Value',0);
                eaApplyCallback(hEdit,hApply,name,FlagNumeric);
            else
                
                % enable the apply button
                set(hApply,'enable','on');
            end
        end
    end % END function eaEditCallback

    function eaEditKeypress(evt,hEdit)
        % EAEDITKEYPRESS Generic edit keypress for a specific context
        %
        %   EA stands for Edit / Apply.  It describes a specific setup of
        %   UI controls in which a user may change an edit box value, and
        %   either press enter or click an apply button.
        %
        %   This keypress's purpose is to detect when an enter or return
        %   key is pressed.
        
        if strcmpi(evt.Key,'enter')||strcmpi(evt.Key,'return')
            set(hEdit,'Value',1);
        end
    end % END function peaEditKeypress

    function eaApplyCallback(hEdit,hApply,name,FlagNumeric)
        % EAAPPLYCALLBACK Generic button callback for a specific context
        %
        %   EA stands for Edit / Apply.  It describes a specific setup of
        %   UI controls in which a user may change an edit box value, and
        %   either press enter or click an apply button.
        %
        %   This callback's purpose is to save the value once the user has
        %   pressed the button.
        
        % make sure apply button is off
        set(hApply,'Enable','off');
        
        % get the new value
        val = get(hEdit,'String');
        if FlagNumeric
            if strcmpi(val(1),'[')
                val = eval(val);
            else
                val = str2double(val);
            end
        end
        
        % udpate the requested value and restore others to original values
        this.(name) = val;
        updateGUI_DisplayFile(this);
    end % END function eaApplyCallback
end % END function layout