function getRerefFigs(subj,noreref,commonAv,figSaveDir,ecog_sessName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

reref = noreref - commonAv;
displaySec = 2; samplerate = 30000;
% if time series is too short
if length(commonAv) > displaySec*samplerate
    
    fontSize = 14;
    
    % Get rid of underscore
    ecog_sessName(strfind(ecog_sessName,'_')) = '-';
    
    % Plot!
    %randomShift = floor(rand*100000);
    timeRange = floor((length(reref)/2):(length(reref)/2)  +  samplerate*displaySec-1);  % 2 seconds
    %timeRange = timeRange + randomShift;
    % quantileRange = [0.025 0.25 0.50 0.75 0.975];
    quantileRange = [0.025 0.975];
    colors = {'r','b','c'};
    toPlot = {commonAv,noreref,reref};
    labels = {'commonAv','noreref','reref'};
    axesA = []; axesB = []; axesC = []; axesD = []; axesE = [];
    voltageLim = 0; logCountLim = 0;
    masterFig = figure('units','normalized','outerposition',[0,0,1,1],'visible','off');
    for pp = 1:length(toPlot)
        figureHandle = subplot(4,3,pp+(2*pp-2)); hold on;
        title(sprintf('%s',labels{pp}));
        plot(1:samplerate*displaySec,toPlot{pp}(timeRange));
        if voltageLim < max(abs(get(gca,'YLim'))), voltageLim = max(abs(get(gca,'YLim'))); end
        set(gca,'XTick',0:samplerate:samplerate*displaySec,'XTickLabel',0:1:displaySec);
        xlabel('Seconds'); ylabel('Voltage(uV)');
        axesA(pp) = gca; setFigFontSize(fontSize,figureHandle);
    end
    for pp = 1:length(toPlot)
        % Wider version of the distribution plot
        figureHandle = subplot(4,3,pp+(2*pp-1)); hold on;
        [oc,bins] = hist(toPlot{pp},-voltageLim*3/2:2:voltageLim*3/2);
        plot(bins,oc,'Color',colors{pp},'LineWidth',2); set(gca,'yscale','log');
        qVals = quantile(toPlot{pp},quantileRange);
        for q = 1:length(quantileRange)
            plot([qVals(q),qVals(q)],get(gca,'YLim'),'k--','LineWidth',1);
        end
        quantWidth = floor(range(qVals));
        text(qVals(1)+20,median(oc),sprintf('Width(95%%): %d',quantWidth));
        if logCountLim < max(abs(get(gca,'YLim'))), logCountLim = max(abs(get(gca,'YLim'))); end
        xlabel('Voltage(uV)'); ylabel('Log(count)');
        axesB(pp) = gca; setFigFontSize(fontSize,figureHandle);
        if pp==1, title(sprintf('%s-%s',subj,ecog_sessName)); end
        if pp==3, textLoc = [qVals(1)+20,median(oc)]; end
        if pp==2, norerefWidth = quantWidth;
        elseif pp==3, rerefWidth = quantWidth;
        end
        
        % Plot with all three distributions
        figureHandle = subplot(4,3,11); hold on;
        plot(bins,oc,'Color',colors{pp},'LineWidth',2); set(gca,'yscale','log');
        xlabel('Voltage(uV)'); ylabel('Log(count)');
        if pp==length(toPlot)
            text(textLoc(1),textLoc(2),sprintf('W_{reref} / W_{noreref} = %0.2f',rerefWidth/norerefWidth));
        end
        axesD(pp) = gca; setFigFontSize(fontSize,figureHandle);
        
        % Narrower version of distribution plot
        figureHandle = subplot(4,3,pp+(2*pp)); hold on;
        [oc,bins] = hist(toPlot{pp},-voltageLim*5:1:voltageLim*5);
        plot(bins,oc,'Color',colors{pp},'LineWidth',2); set(gca,'yscale','log');
        qVals = quantile(toPlot{pp},quantileRange);
        for q = 1:length(quantileRange)
            plot([qVals(q),qVals(q)],get(gca,'YLim'),'k--','LineWidth',1);
        end
        xlabel('Voltage(uV)'); ylabel('Log(count)');
        axesC(pp) = gca; setFigFontSize(fontSize,figureHandle);
        
        % Plot with all three distributions
        figureHandle = subplot(4,3,12); hold on;
        plot(bins,oc,'Color',colors{pp},'LineWidth',2); set(gca,'yscale','log');
        xlabel('Voltage(uV)'); ylabel('Log(count)');
        axesE(pp) = gca; setFigFontSize(fontSize,figureHandle);
    end
    for pp = 1:length(toPlot)
        set(axesA(pp),'YLim',[-voltageLim,voltageLim]);
        set(axesB(pp),'XLim',[-voltageLim,voltageLim],'YLim',[0,logCountLim]);
        set(axesC(pp),'XLim',[-voltageLim*5,voltageLim*5],'YLim',[0,logCountLim]);
        set(axesD(pp),'XLim',[-voltageLim,voltageLim],'YLim',[0,logCountLim]);
        set(axesE(pp),'XLim',[-voltageLim*5,voltageLim*5],'YLim',[0,logCountLim]);
    end
    masterFig = tightfig(masterFig);
    set(masterFig, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
    
    % Save the figure
    print(masterFig,figSaveDir,'-dpng','-r600');
    
    %- make sure everything is done and then close the fig.
    drawnow();
    pause(2);
    close(masterFig);
    
end

end


function setFigFontSize(fontSize,figureHandle)
% Sets the font size of all text in a currently open figure
set(gca,'FontSize',fontSize);
%     figureHandle = gcf;
set(findall(figureHandle,'type','text'),'FontSize',fontSize);
end


function hfig = tightfig(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
%
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
%
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

if nargin == 0
    hfig = gcf;
end

% There can be an issue with tightfig when the user has been modifying
% the contnts manually, the code below is an attempt to resolve this,
% but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
set(hfig, 'WindowStyle', 'normal');

% 1 point is 0.3528 mm for future use

% get all the axes handles note this will also fetch legends and
% colorbars as well
hax = findall(hfig, 'type', 'axes');

% get the original axes units, so we can change and reset these again
% later
origaxunits = get(hax, 'Units');

% change the axes units to cm
set(hax, 'Units', 'centimeters');

% get various position parameters of the axes
if numel(hax) > 1
    %         fsize = cell2mat(get(hax, 'FontSize'));
    ti = cell2mat(get(hax,'TightInset'));
    pos = cell2mat(get(hax, 'Position'));
else
    %         fsize = get(hax, 'FontSize');
    ti = get(hax,'TightInset');
    pos = get(hax, 'Position');
end

% ensure very tiny border so outer box always appears
ti(ti < 0.1) = 0.15;

% we will check if any 3d axes are zoomed, to do this we will check if
% they are not being viewed in any of the 2d directions
views2d = [0,90; 0,0; 90,0];

for i = 1:numel(hax)
    
    set(hax(i), 'LooseInset', ti(i,:));
    %         set(hax(i), 'LooseInset', [0,0,0,0]);
    
    % get the current viewing angle of the axes
    [az,el] = view(hax(i));
    
    % determine if the axes are zoomed
    iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
    
    % test if we are viewing in 2d mode or a 3d view
    is2d = all(bsxfun(@eq, [az,el], views2d), 2);
    
    if iszoomed && ~any(is2d)
        error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.')
    end
    
end

% we will move all the axes down and to the left by the amount
% necessary to just show the bottom and leftmost axes and labels etc.
moveleft = min(pos(:,1) - ti(:,1));

movedown = min(pos(:,2) - ti(:,2));

% we will also alter the height and width of the figure to just
% encompass the topmost and rightmost axes and lables
figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);

figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);

% move all the axes
for i = 1:numel(hax)
    
    set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
    
end

origfigunits = get(hfig, 'Units');

set(hfig, 'Units', 'centimeters');

% change the size of the figure
figpos = get(hfig, 'Position');

set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);

% change the size of the paper
set(hfig, 'PaperUnits','centimeters');
set(hfig, 'PaperSize', [figwidth, figheight]);
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition',[0 0 figwidth figheight]);

% reset to original units for axes and figure
if ~iscell(origaxunits)
    origaxunits = {origaxunits};
end

for i = 1:numel(hax)
    set(hax(i), 'Units', origaxunits{i});
end

set(hfig, 'Units', origfigunits);

%      set(hfig, 'WindowStyle', origwindowstyle);

end










