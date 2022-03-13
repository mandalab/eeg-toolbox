function utahMap = convertUtahElecs2Chans(utahMap, SHOW_PLOT)
%  function convertUtahElecs2Chans
%
%    takes the subject specific "utahMap" fields, which maps "physical electrode numbers" to "space"
%    and converts this to a map between "recorded channel numbers" and "space", which is required to interpret
%    the spatial relationship of recorded spikes and/or LFPs on the array
%
%  This is a helper function is called by "mapUtahChan2Elec", which populates the user-specific utahMap fields
%
%- Input utahMap data -
% utahMap.subj               = subj;
% utahMap.elecNumSpace1      = elecNumSpace1;
% utahMap.elecNumSpace1_wire = elecNumSpace1_wire;
% utahMap.elecNumSpace2      = elecNumSpace2;
% utahMap.elecNumSpace2_wire = elecNumSpace2_wire;
% utahMap.stimPigtail2Chan   = stimPigtail2Chan;
%
%- Output utahMap struct has these additional fields
% utahMap.chanNumSpace1      = chanNumSpace1;
% utahMap.chanNumSpace2      = chanNumSpace2;
% utahMap.chanNumBPlist      - 2xN list of adjacent channel nummbers for bipolar referencing
% utahMap.stimPigtail2Elec   = stimPigtail2Elec;
%
% JW 2/2019
%


if nargin<2, SHOW_PLOT=1; end;


%- copy columns "Ch#" and "Pin #" from Bank A,B,C,D, incrementing by 32 for each bank
chanNumMap = [
    1	2	33	34	65	66	97	98
    3	4	35	36	67	68	99	100
    5	6	37	38	69	70	101	102
    7	8	39	40	71	72	103	104
    9	10	41	42	73	74	105	106
    11	12	43	44	75	76	107	108
    13	14	45	46	77	78	109	110
    15	16	47	48	79	80	111	112
    17	18	49	50	81	82	113	114
    19	20	51	52	83	84	115	116
    21	22	53	54	85	86	117	118
    23	24	55	56	87	88	119	120
    25	26	57	58	89	90	121	122
    27	28	59	60	91	92	123	124
    29	30	61	62	93	94	125	126
    31	32	63	64	95	96	127	128];


%- do the mapping
chanNumSpace1 = nan(size(utahMap.elecNumSpace1));
chanNumSpace2 = nan(size(utahMap.elecNumSpace2));
chanSpaceInfo = nan(numel(utahMap.elecNumMap),4); %- chanNum, device, x, y
for iC=1:numel(utahMap.elecNumMap),
    thisChan = chanNumMap(iC);
    thisElec = utahMap.elecNumMap(iC);
    iE = find(utahMap.elecNumSpace1==thisElec);
    if length(iE)==1,
        chanNumSpace1(iE) = thisChan;
        [y, x] = ind2sub(size(chanNumSpace1),iE);
        chanSpaceInfo(thisChan,:) = [thisChan 1 x y];
    else
        iE = find(utahMap.elecNumSpace2==thisElec);
        chanNumSpace2(iE) = thisChan;
        [y, x] = ind2sub(size(chanNumSpace2),iE);
        chanSpaceInfo(thisChan,:) = [thisChan 2 x y];
    end
end
stimPigtail2Chan = nan(size(utahMap.stimPigtail2Elec));
for iC = 1:size(utahMap.stimPigtail2Elec,1),
    if ~isnan(utahMap.stimPigtail2Elec(iC,1))
        stimPigtail2Chan(iC,1) = utahMap.stimPigtail2Elec(iC,1); %-pigtail ring number
        stimPigtail2Chan(iC,2) = chanNumMap(find(utahMap.elecNumMap==utahMap.stimPigtail2Elec(iC,2))); %-pigtail ring number
    end
end

%- create an output table with channel names conveying physical space
[~,iSort]     = sort(chanSpaceInfo(:,1));  %- seems like this is already ordered, but sort here just in case.
chanSpaceInfo = chanSpaceInfo(iSort,:);
chanSpaceTable = array2table(chanSpaceInfo,'VariableNames',{'chanNumRec' 'arrayNum' 'PosX' 'PosY'});
for iC=1:height(chanSpaceTable),
    if isempty(chanNumSpace2), devStr = '';
    else                       devStr = sprintf('%da',chanSpaceTable.arrayNum(iC)); end
    if max(chanSpaceTable.PosX)<10,
        chanSpaceTable.PosStr{iC}    = sprintf(  '%dx%dy',       chanSpaceTable.PosX(iC),chanSpaceTable.PosY(iC));
        chanSpaceTable.PosStrDev{iC} = sprintf('%s%dx%dy',devStr,chanSpaceTable.PosX(iC),chanSpaceTable.PosY(iC));
    else
        chanSpaceTable.PosStr{iC}    = sprintf(  '%02dx%02dy',       chanSpaceTable.PosX(iC),chanSpaceTable.PosY(iC));
        chanSpaceTable.PosStrDev{iC} = sprintf('%s%02dx%02dy',devStr,chanSpaceTable.PosX(iC),chanSpaceTable.PosY(iC));
    end
end
     

%- BIPOLAR CHANNELS: create list of adjacent channels by combing across chanNumSpace
bpList = [];
for whichAr=[1 2],
    if whichAr==1,chanNumSpace = chanNumSpace1;
    else          chanNumSpace = chanNumSpace2;
    end
    for iRow=1:size(chanNumSpace,1),
        for iCol=1:size(chanNumSpace,2),
            if iRow<size(chanNumSpace,1), bpList = [bpList; chanNumSpace(iRow,iCol) chanNumSpace(iRow+1,iCol)]; end %- left to right
            if iCol<size(chanNumSpace,2)  bpList = [bpList; chanNumSpace(iRow,iCol) chanNumSpace(iRow,iCol+1)]; end %- top to bottom
        end
    end
end %- for whichAr
notNan = find(~isnan(sum(bpList,2))); %- 96 channel arrays intentionally have nans in the corners... cut those out of the bp list
bpList = bpList(notNan,:);


%- organize the list.  first make it so smaller value is always in first column, the ascending 2nd column, then ascending 1st column

bpListSort = sort(bpList')';       %- lower index in fist column
[~,iSort]  = sort(bpListSort(:,2)); %- order the second column
bpListSort = bpListSort(iSort,:);
[~,iSort]  = sort(bpListSort(:,1)); %- order the first column
bpListSort = bpListSort(iSort,:);
bpList     = bpListSort;





%- some nans found... why?
%x = sort(elecNumMap(:));
%y = diff(x);
%[x [y; 1]];


%infoStr = sprintf('%s\n subj               = subject string (e.g., NIH029)',utahMap.infoStr);
infoStr = utahMap.readme;
infoStr = sprintf('%s\n chanNumSpace1      = maxtrix (num utah rows x col), of recorded channel numbers mapped to physical space',infoStr);
infoStr = sprintf('%s\n chanNumSpace2      = if multi-port (two arrays), maxtrix (num utah rows x col) for second array',infoStr);
infoStr = sprintf('%s\n chanSpaceTable     = table (# electrodes x 6), with fields: chanNumRec arrayNum PosX PosY PosStr PosStrDev ',infoStr);
infoStr = sprintf('%s\n           where: PosX=1 and PosY=1 refered to the upper left corner of the chanNumSpace matrix. ',infoStr);
infoStr = sprintf('%s\n           where: PosStr contains strings intended to be used in channel names for the array ("XxYy"), PosStrDev additionally contains array number',infoStr);
infoStr = sprintf('%s\n bpList             = matrix (# bp pairs x 2), of all adjacent micro electrodes, where each value is an recorded channel num',infoStr);
infoStr = sprintf('%s\n stimPigtail2Chan   = map from stimulation channels (1-15) to recorded channle numbers',infoStr);


%- save the mapped data to the utah structure
utahMap.readme           = infoStr;
utahMap.chanNumSpace1    = chanNumSpace1;
utahMap.chanNumSpace2    = chanNumSpace2;
utahMap.chanSpaceTable   = chanSpaceTable;
utahMap.chanNumBPlist    = bpList; %- list of channel numbers corresponding to adjacent electrodes.
utahMap.stimPigtail2Chan = stimPigtail2Chan;






%- Make a graphic showing where everything is
if SHOW_PLOT,
    
    %- multi-port (2 arrays) or single array
    if length(utahMap.elecNumSpace2)==0,
        multiPort=0; iPlotSet=[1 3]; figWidth = 550;
    else
        multiPort=1; iPlotSet=[1:4]; figWidth = 1100;
    end
    
    
    figure(1); clf
    set(gcf,'color','w','name',sprintf('%s Utah Electrode Mapping',utahMap.subj));
    set(gcf,'position', [100         350        figWidth         950]);
    fS = 15;
    
    
    for iPlot=iPlotSet,
        
        eSpaceWire=[];
        switch iPlot,
            case 1,
                eSpace     = utahMap.elecNumSpace1;
                eSpaceWire = utahMap.elecNumSpace1_wire;
                stimMap    = utahMap.stimPigtail2Elec(:,2);
                tStr = 'Array 1: Elec Nums in space';
            case 2
                eSpace     = utahMap.elecNumSpace2;
                eSpaceWire = utahMap.elecNumSpace2_wire;
                stimMap    = utahMap.stimPigtail2Elec(:,2);
                tStr = 'Array 2: Elec Nums in space';
                
            case 3
                eSpace     = utahMap.chanNumSpace1;
                for iC=1:length(utahMap.elecNumSpace1_wire), eSpaceWire = [eSpaceWire utahMap.chanNumSpace1(find(utahMap.elecNumSpace1(:)==utahMap.elecNumSpace1_wire(iC)))]; end
                stimMap    = utahMap.stimPigtail2Chan(:,2);
                tStr = 'Array 1: Chan Nums in space';
                
            case 4
                eSpace     = utahMap.chanNumSpace2;
                for iC=1:length(utahMap.elecNumSpace2_wire), eSpaceWire = [eSpaceWire utahMap.chanNumSpace2(find(utahMap.elecNumSpace2(:)==utahMap.elecNumSpace2_wire(iC)))]; end
                stimMap    = utahMap.stimPigtail2Chan(:,2);
                tStr = 'Array 2: Chan Nums in space';
        end
        
        if multiPort,
            subplot(2,2,find(iPlot==iPlotSet))
        else
            subplot(2,1,find(iPlot==iPlotSet))
        end
        nRow = size(eSpace,1);
        nCol = size(eSpace,2);
        
        imagesc([1:nRow],[1:nCol],eSpace,[1 max([utahMap.chanNumSpace1(:); utahMap.chanNumSpace2(:)])]); hold on;
        set(gca,'fontsize',fS,'xtick',[1:nCol],'xticklabel',{},'ytick',[1:nRow],'yticklabel',{},'box','off','tickdir','out');
        title(tStr);
        colormap(cool)
        axis square
        %-
        for iRow=1:nRow,
            for iCol=1:nCol,
                if     any(eSpace(iCol,iRow)==eSpaceWire),  bg = 'r';
                elseif any(eSpace(iCol,iRow)==stimMap),     bg = 'y';
                else                                        bg = 'w';       end
                hT1 = text(iRow,iCol,sprintf('%d',eSpace(iCol,iRow)),'fontsize',fS,'color','k', 'FontWeight','bold',  'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',bg);
            end
        end
        
        if iPlot==1,
            xlabel('RED = electrodes neighboring wire bundle');
        elseif iPlot==2,
            xlabel('YELLOW = stim electrodes');
        end
        
    end
    hT = suplabel(utahMap.subj  ,'t');
    set(hT,'fontsize',fS+5)
    
end % if SHOW_PLOT
