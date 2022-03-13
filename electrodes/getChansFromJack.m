function [jschans, physicalChan, nums] = getChansFromJack(jacksheet_in)
%getChansFromJack is a function to centralize the pulling of channel names from jacksheetBR/CV_local,
%which happens in split_BR, split_CV, and getChanNames_CV
%
% The user defines either the path to jacksheetBR/CV_local.csv or the same file read in as a table.
% getChansFromJack will determine whether this is Blackrock (BR) or Cervello (CV) data and isolate the
% same channels that would be output in ChanNames.txt (now obsolete) and used for splitting.
%
% INPUT:
%     jacksheet_in (char OR table): the path (char) to a jacksheetBR/CV_local.csv, OR the same file
%                                   already read in with readtableSafe (table)
% 
% OUTPUTS:
%     jschans                     : (cell array) list of channels
%     physicalChan                : (cell array) list of physicalChans (only for CV)
%     nums                        : (double) number of channels (only for BR)
% 
%
%
% Get the chanNames as would be expected in the obsolete ChanNames.txt, from jacksheetBR_local or
% jacksheetCV_local, OR a MATLAB variable
% 
% Called by getChanNames_CV (not BR version b/c NSP-dependent)
%
% Created by Samantha Nicole Jackson, 9/29/2020
%
% REVISION HISTORY:
% 09/29/2020 SJ: created function
% 10/28/2020 SJ: added din to list of what to exclude from BR (shows up in nev cases?)
% 01/22/2021 SJ: added ain21 to exclusion list 
%
% To DO: If the jacksheet is _complete, then verify that is being called by phys2name, and maybe we will
% have to do something different (because no chanNameNew column!!)

% First check to see if the input is a char or a table
if ischar(jacksheet_in)
    if contains(jacksheet_in,'jacksheetBR_local.')
        hwtype = 'BR';
    elseif contains(jacksheet_in,'jacksheetCV_local.')
        hwtype = 'CV';
    else
        fprintf('%s\n',['ERROR!!! jacksheetPath does not contain jacksheetBR_local. or jacksheetCV_local. (' jacksheet_in ') Get SJ!']);
        keyboard
    end
    jacksheet = readtableSafe(jacksheet_in);
    
elseif istable(jacksheet_in)
    % Figure out if it is BR or CV
    if any(strcmpi('NSPsuffix',jacksheet_in.Properties.VariableNames)) %Only a BR field
        if size(jacksheet_in,2) <= 4
            keyboard %Only CV jacksheets should have 4 or less columns! Yet it has NSPsuffix?
        end
        hwtype = 'BR';
    else
        if size(jacksheet_in,2) > 4
            fprintf('%s\n','ERROR!! Looks like jacksheetBR_complete is being used because jacksheetCV_local has not been created!');
            keyboard %Only BR jacksheets should have more than 4 columns! Yet it does not has NSPsuffix?
        end
        hwtype = 'CV';
    end
    jacksheet = jacksheet_in;
else
    fprintf('%s\n','ERROR!! Input jacksheet_in is not a char path to a jacksheet or a jacksheet table! Get SJ');
    keyboard
end

% Get the channels
if strcmpi(hwtype,'BR')
    jschans = cell(sum(~contains(jacksheet.FileName,'.nev') & cellfun(@isempty,regexp(jacksheet.ChanName,'(ain.*(17|18|19|20|21))|din1.*'))),1); %ain21 added 1/22/2021
    %SJ: get nums now by searching for the 1st NSP but not the pulse channels,
    %which should hopefully be at the end...
    nums = sum(jacksheet.NSP == 1 & jacksheet.PhysicalChan < 1000);
    ii_row = 0;
    for ii = 1:numel(jacksheet.ChanName)
        if ~contains(jacksheet.FileName{ii},'.nev') && isempty(regexp(jacksheet.ChanName{ii},'(ain.*(17|18|19|20|21))|din1.*')) %SJ added to take care of .ns2 ain17-20s; later added din for nev cases..?; later for ain21
            ii_row = ii_row + 1;
            if strcmp(jacksheet.ChanNameNew{ii},'-')
                jschans{ii_row} = jacksheet.ChanName{ii};
            else
                jschans{ii_row} = jacksheet.ChanNameNew{ii};
            end
        end
    end
    physicalChan = {}; % Not for BR data, this is handled in getChanNames_BR
    
elseif strcmpi(hwtype,'CV')
    jschans = cell(numel(jacksheet.ChanName),1);
    physicalChan = cell(numel(jacksheet.ChanName),1);
    for jj = 1:numel(jacksheet.ChanName)
        if strcmp(jacksheet.ChanNameNew{jj},'-')
            jschans{jj} = jacksheet.ChanName{jj};
        else
            jschans{jj} = jacksheet.ChanNameNew{jj};
        end
        physicalChan{jj} = sprintf('PhysicalChan=%d',jacksheet.PhysicalChan_guess(jj));
    end
    % Now we need to remove the 'exclude' channels (din*, ain17-20) -> from
    % getChanNames_CV (even though these are written to jacksheetCV)
    exclude_chans = ~ismember(jschans,regexpCell(jschans,'^(din.*|ain1[7-9]|ain2[0-1])$')); %SJ added ain21
    jschans = jschans(exclude_chans);
    physicalChan = physicalChan(exclude_chans);
    nums = 0; % Not for CV data
else
    keyboard %Not BR or CV??
end


end

