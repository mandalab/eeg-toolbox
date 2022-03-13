function [microNameTable, jackTable] = microNameConvert(subjectInfo,jackTable,SKIP_KEYBOARD)
%- function microNameConvert

%- input:  subject       -- 'NIH029'
%-         tJacksheetBR  -- table from jacksheetBR_complete
%                 code below uses fields: ChanName, NSPsuffix, PhysicalChan, SampFreq, FileName
%                  to update/adds fields: ChanNameNew, MicroDevNum
%
%  output:  tNameMicro   -- table:   NSPsuffix, oldChanStr, newChanStr, devNum
%        jackTableUpdate -- table: original jacktable, but with NSPsuffix, MicroDevNum and ChanNameNew added/updated
%
%
%  created 3/2019 by JHW
%

if nargin<3,
    SKIP_KEYBOARD = 0; %- assume you cannot skip.  Only SKIPs keyboards that were *warnings*... *errors* still cause a keyboard break
end

%%- TEST INPUTS
%subjectInfo    = getMicroSubjInfo_v11('NIH069');
%jackTable      = readtable('/Volumes/JW24TB/data24TB/localFRNU56/NIH069/data_raw/190125_1045/jacksheetBR_complete.csv');


i30k           = find(jackTable.SampFreq==30000 & ~contains(jackTable.FileName,'.nev') & jackTable.PhysicalChan<1000); %- MOSTLY micros... in early subjects (<~NIH040) occasionally used on eCog and/or stim DC signals
jacksheet30kHz = jackTable(i30k,:);

%- make a temporary column tagName... for each ChanName create a string that is missing the numbers
for iRow=1:height(jacksheet30kHz),
    jacksheet30kHz.TagName{iRow} = jacksheet30kHz.ChanName{iRow}(1:find(isstrprop(jacksheet30kHz.ChanName{iRow},'digit'),1,'first')-1);
end

tableCell      = {}; %- initailize a cell array that will get converted to table after filling


%- loop over devices... whether separate utah arrays, microwire depths, etc.  Find the recorded channel that matches "OLD CHAN" from subjectInfo then create a mapping to "NEW CHAN"
nonMatchCount=0;  jacksheetUsed = zeros(height(jacksheet30kHz),1);
for iDEV=1:length(subjectInfo.utahChanCnt)
    
    %- big first check... does this nspSuffix even exist?  might not if the array died and the nsp was turned off
    thisDevSuffix  = subjectInfo.nspSuffix{iDEV};
    %iChan30thisNSP = find( strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix) & contains(jacksheet30kHz.ChanName,subjectInfo.oldChanName{iDEV,1}));
    iChan30thisNSP = find( strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix) & strcmp(jacksheet30kHz.TagName,subjectInfo.oldChanName{iDEV,1})); %- this will help knock out some error; e.g. NIH060 tag LPHD vs LPHDm
    
    %- if not found, see if alt NSP suffix is in play
    if length(iChan30thisNSP)==0 & length(subjectInfo.nspSuffixAlt)>=iDEV,
        thisDevSuffix  = subjectInfo.nspSuffixAlt{iDEV};
        %iChan30thisNSP = find( strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix) & contains(jacksheet30kHz.ChanName,subjectInfo.oldChanName{iDEV,1}));
        iChan30thisNSP = find( strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix) & strcmp(jacksheet30kHz.TagName,subjectInfo.oldChanName{iDEV,1})); %- this will help knock out some error; e.g. NIH060 tag LPHD vs LPHDm
    end
    
    %- channel x suffix not found...
    if length(iChan30thisNSP)==0,
        %- channel x suffix not found.  IF name or suffix got messed up and dont match other devices we will trigger a flag below
        fprintf('\n NSP suffix "%s" with chan name "%s" not found in session "%s" jacksheet (and chanName not found otherwise), so device "%s" assumed not present', thisDevSuffix,subjectInfo.oldChanName{iDEV,1},jackTable.RawDir{1},subjectInfo.newChanName{iDEV});
        nonMatchCount = nonMatchCount+subjectInfo.utahChanCnt(iDEV); %- definitely want to trigger a flag below that there is an issue (only checked if unaccounted 30k channels)
        continue
    end
    
    
    %- number of channels expected in this micro device (96 or 64+64 for utah; usually 16 for microwires)
    thisCnt = subjectInfo.utahChanCnt(iDEV);
    physChanStart = nan; %- for the first channel of a device, store the "phys chan" number from the jacksheet. This is a unique identifier that can be used to guess the correct subsequent channnels if ChanStr doesnt match expectation
    for iC=1:thisCnt
        
        %% Construct "OLD NAME" from the subjectInfo data and search for a single match in jacksheetBR
        % make sure that oldChanName is a match with the electrode names in the jacksheet.
        % Most likely culpret for a mismatch is leading zeros in number part of the string
        oldName = sprintf('%s%03d',subjectInfo.oldChanName{iDEV,1},subjectInfo.oldChanName{iDEV,2}+iC-1);
        firstTry = oldName; %- for potential error message below
        matchWithRef = find(strcmp(jacksheet30kHz.ChanName,oldName) & strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix));
        
        if length(matchWithRef) == 0,
            % try different number of 0s. probably not the most elegant way to do this
            oldName = sprintf('%s%02d',subjectInfo.oldChanName{iDEV,1},subjectInfo.oldChanName{iDEV,2}+iC-1);
            matchWithRef = find(strcmp(jacksheet30kHz.ChanName,oldName)& strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix));
            if length(matchWithRef) == 0,
                oldName = sprintf('%s%01d',subjectInfo.oldChanName{iDEV,1},subjectInfo.oldChanName{iDEV,2}+iC-1);
                matchWithRef = find(strcmp(jacksheet30kHz.ChanName,oldName) & strcmp(jacksheet30kHz.NSPsuffix,thisDevSuffix));
                if length(matchWithRef) == 0,
                    % even with changing 0s, this chan is not in any of the files to be extracted
                    fprintf(sprintf('\nchan name %s [%s] from getMicroSubjInfo doesnt match any channel in jacksheets',oldName,thisDevSuffix));
                    nonMatchCount=nonMatchCount+1;
                end
            end
        end
        
        %- final check... should be exactly 1 match
        if     length(matchWithRef)==1 & jacksheetUsed(matchWithRef)==0, jacksheetUsed(matchWithRef)=1 ;
        elseif length(matchWithRef)==1 & jacksheetUsed(matchWithRef)>0,  fprintf('\n ERROR: shouldnt get here...'); keyboard; jacksheetUsed(matchWithRef)=jacksheetUsed(matchWithRef)+1 ;
        elseif length(matchWithRef)==0, oldName = 'ChanNotFound';
        elseif length(matchWithRef)>1,  fprintf('\n ERROR: shouldnt get here. fix the issue'); keyboard; end
        
        
        
        %% USE THE PHYSICAL CHANNEL NUMBERS TO CHECK/IDENTIFY "OLD CHANNELS"
        %- if first channel of this device record the phys num.  If subsequent channels, confirm found physChanNum matches expectation. Also use this info to fix the channels that dont matchRef
        if iC==1,
            physChanStart = jacksheet30kHz.PhysicalChan(matchWithRef);
        end
        physChanExpected  = physChanStart+iC-1;
        if ~isempty(physChanExpected),
            matchPhysChan = find(jacksheet30kHz.PhysicalChan == physChanExpected);
        else
            matchPhysChan = []; 
        end
            
        if length(matchWithRef)==1,  %- sanity check... matchWithRef SHOULD equal matchPhysChan
            if matchPhysChan~=matchWithRef,
                fprintf('\n ERROR: unexpected mismatch between expected and actual phys chan. figure out why. \n note: code below assumes physChanExpected equals physChan of this channel\n');
                keyboard;
            end
            
        elseif length(matchPhysChan)==1, %- if oldChan name was not found, last possible fix would be to use expected physChan number to pick the channel.  Try it here
            if jacksheetUsed(matchPhysChan)>0,
                fprintf('\n ERROR: channel %s not found.channel %s is in the correct spot in jacksheet, but that spot has already been used. Not sure what to do',firstTry,jacksheet30kHz.ChanName{matchPhysChan});
                keyboard;
            else
                fprintf('\n');
                iShow = find(jacksheet30kHz.PhysicalChan >= physChanExpected-2 & jacksheet30kHz.PhysicalChan <= physChanExpected+2);
                disp(jacksheet30kHz(iShow,:));
                if 1 | SKIP_KEYBOARD,
                    fprintf('\n HEADS UP: channel %s not found, but channel %s is in the correct spot in the jacksheet and will be used in its place', firstTry,jacksheet30kHz.ChanName{matchPhysChan});
                    reply='y';
                else
                    [reply] = input(sprintf('\n WARNING: channel %s not found, but channel %s is in the correct spot in the jacksheet\n and that spot has not been used... use it?',firstTry,jacksheet30kHz.ChanName{matchPhysChan}),'s');
                    if isempty(reply), reply = 'n'; end
                end
                if strcmpi(reply,'y'),
                    matchWithRef = matchPhysChan;
                    oldName     = jacksheet30kHz.ChanName{matchWithRef};
                    jacksheetUsed(matchWithRef)=jacksheetUsed(matchWithRef)+1;
                    nonMatchCount = nonMatchCount-1; %- since we commited to the name change, dont trigger a warning below about this channel
                end
                
            end
            
        else %- oldChanName not found, and physChan isn't pointing to the channel... looks like it wasn't recorded
            fprintf('\n ERROR: channel "%s" not found at 30kHz. New name would have been "%s". Appears to not have been recorded. ',firstTry,subjectInfo.newChanName{iDEV,1});
            %- this error is triggered with NIH060 because device exclusion above is trigged by a very similar name (LPHDm vs LPHD... contains hits for both)
            if sum(jacksheetUsed)<length(jacksheetUsed)
                fprintf('\n     and there are unused channels in jacksheet at this point, so worth taking a closer look');
                keyboard
            end
        end
        
        
        %% ANOTHER DATA CHECK.  MAKE SURE actual NSP matches subjectInfo.nspSuffix.  Above only checked in case of two channel names
        thisNSPsuffix = '';
        if length(matchWithRef==1),
            if ~contains(jacksheet30kHz.FileName{matchWithRef}, thisDevSuffix),
                fprintf('\n WARNING: jacksheet nspFileName (%s) does not match NSPsuffix (%s) from getMicroSubjInfo... check into this',jacksheet30kHz.FileName{matchWithRef}, thisDevSuffix);
                keyboard;
            elseif ~strcmp(jacksheet30kHz.NSPsuffix{matchWithRef}, thisDevSuffix),
                %- slightly more stringent requirement/check... that NSPsuffix matches exactly
                fprintf('\n WARNING: jacksheet NSPsuffix (%s) does not match NSPsuffix (%s) from getMicroSubjInfo... check into this',jacksheet30kHz.NSPsuffix{matchWithRef}, thisDevSuffix);
                keyboard;
            else
                thisNSPsuffix = jacksheet30kHz.NSPsuffix{matchWithRef};
            end
        end
        
        
        %- Now create the "new channel name".   This comes from getMicroSubjInfo, which calls getUtahChannelMap to get physical XY locations for each channel.
        if contains(subjectInfo.newChanName{iDEV,1},'utah'),
            utahRecChan = mod(physChanExpected-1,128)+1;
            utahLocStr  = subjectInfo.utahMap.chanSpaceTable.PosStr{utahRecChan};
            %utahLocStr  = subjectInfo.utahMap.chanSpaceTable.PosStrDev{utahRecChan}; %- longer name with array number
            
            %newName = sprintf('%s_%s',subjectInfo.newChanName{iDEV,1},utahLocStr);   %- channel name based on electrode LOCATION
            newName = sprintf('%s[%s]',subjectInfo.newChanName{iDEV,1},utahLocStr);   %- channel name based on electrode LOCATION
            
            
            %- Old way of doing it.. based on recorded channel number
            %newName = sprintf('%s%02d',subjectInfo.newChanName{iDEV,1},subjectInfo.newChanName{iDEV,2}+iC-1);  %- name it based on the recording channel number
        else
            %- for depths, just rename to match the recording order.  One day this could get converted to linear x radial notation for microwire depths
            newName = sprintf('%s%02d',subjectInfo.newChanName{iDEV,1},subjectInfo.newChanName{iDEV,2}+iC-1);
        end
        
        %- make the sort output name linked to the NSPsuffix and chan name, that way if a new naming scheme comes out the sorts dont have to be redone
        sortChanName =  sprintf('[%s]%s',thisNSPsuffix,oldName);
        
        iT = size(tableCell,1)+1;
        tableCell{iT,1} = thisNSPsuffix;
        tableCell{iT,2} = oldName;
        tableCell{iT,3} = sortChanName;
        tableCell{iT,4} = newName;
        tableCell{iT,5} = subjectInfo.deviceNum(iDEV);
        
    end
end

if isempty(tableCell),
%     if height(jacksheet30kHz)>0, %- will get checked on in a moment
%         fprintf('\n');
%         disp(jacksheet30kHz)
%         fprintf('\n Got to the end of microName without any matches... Does this really happen?');
%         if SKIP_KEYBOARD,  fprintf('\n pausing for 5 sec'); pause(5); else keyboard; end
%     end
    tableCell = {'','','','',[]};
end

microNameTable = cell2table(tableCell,'VariableNames',{'NSPsuffix' 'ChanName' 'SortChanName' 'ChanNameNew' 'MicroDevNum'});



%%% SANITY CHECKS
%  #1:  ARE THERE 30k channels recorded that did not match up with micro?  List set of exclusions for known subjects/channels falling into this category
if sum(jacksheetUsed)<length(jacksheetUsed),
    
    
    %%  list of exclusions for each subject... examples of chan name, suffix, and subject that were recorded at 30k but we know should not be considered micro channels
    %- assume unused 30k channels this is an error that should be reported... unless it fits into one of the exclusions listed below and can be skipped
    
    iExclude = [];
    
    switch subjectInfo.subj,
        
        case 'NIH029'  %-
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.ChanName,'ieeg') & strcmp(jacksheet30kHz.NSPsuffix,'ieeg'));
            if length(iExclude)>0, fprintf('\n HEADS UP: unused 30k channels ieeg from NIH029 appear to be eCog at high sample rate.'); end

        case 'NIH039'  %- 
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.ChanName,'ieeg') & strcmp(jacksheet30kHz.NSPsuffix,'utah_m'));
            if length(iExclude)>0, fprintf('\n HEADS UP: unused 30k channels ieeg from NIH039 appear to be eCog or Noise after array failed at high sample rate.'); end

        case 'NIH042'  %- 
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.ChanName,'ieeg') & strcmp(jacksheet30kHz.NSPsuffix,'utah_m'));
            if length(iExclude)>0, fprintf('\n HEADS UP: unused 30k channels from NIH042 appear to be eCog at high sample rate.'); end
            
        case 'NIH047'
            iExclude = find( jacksheetUsed==0 & (contains(jacksheet30kHz.ChanName,'chan') | strncmp(jacksheet30kHz.ChanName,'el',2)) & strcmp(jacksheet30kHz.NSPsuffix,'ieeg'));
            if length(iExclude)>0, fprintf('\n HEADS UP: unused 30k channels from NIH047 appear to be eCog at high sample rate.'); end
            
        case 'NIH049'
            iExclude = find( jacksheetUsed==0 & strcmp(jacksheet30kHz.ChanName,'LAD4') & strcmp(jacksheet30kHz.NSPsuffix,'micro'));
            if length(iExclude)>0,  fprintf('\n HEADS UP: unused 30k channel from NIH049 (LAD4) was macro contact sampled at 30k for referencing.'); end
            
        case 'NIH053'
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.NSPsuffix,'INST0') );
            if length(iExclude)>0,   fprintf('\n HEADS UP: unused 30k channels from NIH053 appear to be eCog at high sample rate.'); end
   
        case 'NIH062'
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.NSPsuffix,'INST1') & jacksheet30kHz.NSPchan<129);
            if length(iExclude)>0,   fprintf('\n HEADS UP: unused 30k channels from NIH062 appear to be eCog on 2nd NSP at high sample rate.'); end
    
        case 'NIH063'
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.RawDir,'180921_1536'));
            if length(iExclude)>0,   fprintf('\n HEADS UP: unused 30k channels from NIH063 were probably utah recorde with eCog names. only happens once so let it slide.'); end

        case 'NIH071'
            iExclude = find( jacksheetUsed==0 & contains(jacksheet30kHz.ChanName,'microOff') & contains(jacksheet30kHz.RawDir,'190410_1411'));
            if length(iExclude)>0,   fprintf('\n HEADS UP: unused 30k channels from NIH071 were named "microOff" in config as a placeholder when but then were erroneously recorded.'); end
            keyboard %- this happend for one session... any others?
            
    end
    
    
    jacksheetUsed(iExclude)=1; %- apply the exclusion to the jacksheetUsedList
    
    
    
    
%     if sum(jacksheetUsed)<length(jacksheetUsed),
%         fprintf('\n Channels recorded at 30kHz were not specified in getMicroSubjInfo.  Rare but does sometimes happen. Make sure not a mistake:\n');
%         disp(jacksheet30kHz(~jacksheetUsed,:));
%         if SKIP_KEYBOARD,  fprintf('\n pausing for 5 sec'); pause(5); else keyboard; end
%     end
    
    
    %%% SANITY CHECK #2:  ARE THERE MICRO CHANNELS SPECIFIED in MICRO_SUBJ_INFO THAT WERE NOT RECORDED?
    if nonMatchCount>0 & sum(jacksheetUsed)<length(jacksheetUsed),
        %- if you get here, then (1) we didnt find an expected micro device and (2) there are unaccounted channels. DEFINITELY need to confirm that the channels werent just mis-labeled for this session
        fprintf('\n\n Channels specified in getMicroSubjInfo were never found in jacksheets (see above errors), though there are unclaimed 30kHz channels. \nMight need to add another chanNameOld to microSubj. Check into this. Ask JW or Mo\n');
        disp(jacksheet30kHz)
        disp(jacksheet30kHz(~jacksheetUsed,:));
        if SKIP_KEYBOARD, fprintf('\n pausing for 5 sec'); pause(5); else keyboard; end
    end
    
    
elseif nonMatchCount>0,
    %- acceptable to get here.
    fprintf('\n Channels or device specified in getMicroSubjInfo not found, but no unclaimed 30k channels, so probably device was turned off and everything seems to be in balance');
    
    
end






%- now package the result back into the original jackTable
for iRow = 1:height(jackTable),
    iMicro = find( strcmp(microNameTable.ChanName, jackTable.ChanName{iRow}) & strcmp(microNameTable.NSPsuffix, jackTable.NSPsuffix{iRow}) );
    if length(iMicro)==1,
        jackTable.ChanNameNew{iRow} = microNameTable.ChanNameNew{iMicro};  %- do we want to apply the micro name to 1kS data?  I think so...
        jackTable.MicroDevNum(iRow) = microNameTable.MicroDevNum(iMicro);
    else
        jackTable.MicroDevNum(iRow) = -1; %- cant mix strings and numbers... seems better to just make this zero or neg 1
    end
end

