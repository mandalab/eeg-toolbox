function nknih_anon(EEG_file)
% This function anonymizes a Nihon Kohden .EEG file.  It does three things removes the patient's name from the EEG file
% REVISION HISTORY
% 08/31/16 MST - change input parameter from directory to eeg full filename
% 03/09/18 MST - Skip non-EEG files (e.g. blackrock)

[~,~,ext] = fileparts(EEG_file);
if ~(strcmpi(ext, 'EEG') | strcmpi(ext, '.EEG'))  %- jw found that parts is returning .EEG as of 12/2018... how did this ever work?
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .EEG: print and write over name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(exist(EEG_file, 'file') > 0);
fid     = fopen(EEG_file,'r+');

%-read patient ID and overwrite with XXX's
fseek(fid,48,'bof'); 
patID  = fread(fid,16,'*char')';
%XXid  = repmat('X',1,16); 
XXid   = [' ' repmat('X',1,14) ' ']; %add space at begining and end of string
fseek(fid,48,'bof');
fwrite(fid,XXid,'char'); %overwrite the patient ID in the device block

%-read patient name and overwrite with XXX's
fseek(fid,79,'bof'); 
patName = fread(fid,32,'*char')';
%XXname = repmat('X',1,32); 
XXname  = [' ' repmat('X',1,30) ' ']; %add space at begining and end of string
fseek(fid,79,'bof');
fwrite(fid,XXname,'char'); %overwrite the name in the device block
fclose(fid);

if ~strcmp(patID,XXid) || ~strcmp(patName,XXname)
    fprintf('  anonymizing raw .EEG file: removed NAME=%s and ID=%s.\n',patName,patID);
else
    fprintf('  anonymizing raw .EEG file: file already clean.\n');
end
    

%%%%%%%%%%%%%%%%%%
%Completion status
%%%%%%%%%%%%%%%%%%
OUTPUT_STATUS = 0;
if OUTPUT_STATUS
    fprintf('\nDone anonymizing. Run again to make sure names are gone.\n')
    fprintf('OK to copy %s and its .21E to server.\n\n',EEG_file);
end