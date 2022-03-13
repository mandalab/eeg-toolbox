function searchAndDestroyLockedFiles(rootEEGdir)
%
%
%  unlock files on mac os
%
%

%- search for locked files
%cmdSearch = 'find /Volumes/JW24TB/data24TB/eeg_new/ -type f -flags +uchg';
%cmdSearch = 'find /Volumes/jwGDRIVE/data24TB/eeg_new/ -type f -flags +uchg';
%cmdSearch = 'find /Volumes/JW24TB/dataEEGpristine/update/eeg/ -type f -flags +uchg';
%cmdSearch = 'find /Volumes/JW24TB/dataEEGpristine/pristine/eeg/ -type f -flags +uchg';
%cmdSearch = 'find /Volumes/EEG_72PRISTINE/pristine/eeg/ -type f -flags +uchg';
cmdSearch = 'find /Volumes/Shares/FRNU/data/eeg/ -type f -flags +uchg';

cmdSearch = sprintf('find %s -type f -flags +uchg',rootEEGdir);

fprintf('\n\n Running search command: %s',cmdSearch);
[status,result] = system(cmdSearch); 
if length(result)==0, fprintf('\n no locked files found in %s\n',cmdSearch); return; else disp(result); end

%- parse the result
strList = {};
iStr = strfind(result,'/Volumes');
for ii=1:length(iStr)-1,
    thisStr = result(iStr(ii):iStr(ii+1)-2);
    fprintf('\n%d) %s', ii,thisStr);
    strList{end+1} = thisStr;
end
thisStr = result(iStr(end):end-1);
fprintf('\n%d) %s', ii,thisStr);
strList{end+1} = thisStr;

fprintf('\n paste the following line into the terminal so you can respond with password\n\%s\nthen F5 to continue',sprintf('sudo chflags nouchg %s',strList{1}));
system('open ~/../../Applications/Utilities/Terminal.app');
keyboard;

%- unlock each one
for ii=1:length(strList),
    cmd = sprintf('sudo chflags nouchg %s',strList{ii});
    fprintf('\n%d) attempt unlock --> %s', ii, cmd);
    [status,result] = system(cmd);
    fprintf(' status = %d',status);
    %- if above command hangs, then copy first cmd line to terminal,
    %- which will prompt for password.  Once that is in, then rerun loop below and it shoudl be fine
end


[status,result] = system(cmdSearch); 
if length(result)==0, fprintf('\n no locked files found in %s',cmdSearch); return; else disp(result); end


return;


%% List manually copied over from sync error log
strList = {};
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_2/.1032017_06_15_14_39.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_2/.0822017_06_15_14_39.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_3/2017_06_15_14_55.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_3/.1032017_06_15_14_55.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_0/.smbdeleteAAA150903fc9';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_0/2017_06_15_09_27.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_0/.1042017_06_15_09_27.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_1/2017_06_15_11_23.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_1/.1032017_06_15_11_23.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/session_1/.1032017_06_15_11_23.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/junk/session_0/.0822017_06_15_09_25.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/junk/session_jk/2017_06_15_14_37.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH049/behavioral/stimMapGUI/oldVersion/junk/session_jk/.1022017_06_15_14_37.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/session_6/2017_07_12_17_43.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/session_6/.1852017_07_12_17_43.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_3/2017_07_12_17_33.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_3/.1012017_07_12_17_33.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_0/2017_07_12_13_17.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_0/.1042017_07_12_13_17.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_4/2017_07_12_17_38.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_4/.1042017_07_12_17_38.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_5/.smbdeleteAAA150903fd0';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_5/2017_07_12_17_40.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_5/.1032017_07_12_17_40.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_5/.1032017_07_12_17_40.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_2/2017_07_12_17_29.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_2/.1822017_07_12_17_29.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_2/.1822017_07_12_17_29.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/session_7/2017_07_12_17_45.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_2/.1822017_07_12_17_29.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/session_1/2017_07_12_13_23.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/_junk/session_2/.1822017_07_12_17_29.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH050/behavioral/stimMapGUI/session_8/2017_07_12_18_01.pulseLog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH051/behavioral/stimMapGUI/_cantalign_session_10/.smbdeleteAAA1506f434a';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH051/behavioral/stimMapGUI/_cantalign_session_10/2017_07_17_13_54.stimlog';
strList{end+1} = '/Volumes/Shares/FRNU/data/eeg/NIH051/behavioral/stimMapGUI/_cantalign_session_10/.1002017_07_17_13_53.pulseLog';


