function [outputFileName] = bankrelabeling(fileToRead)

%NS5 = openNSx(fileToRead);

%- use concatOpenNSx so there is no cell array of data... that will also make it so a single Header entry is expected (intead of a vector)
NS5 = concatOpenNSx(fileToRead,0,1,[]);
if isempty(NS5.RawData.DataHeader),
    NS5.RawData.DataHeader = 999; %- no idea what should be here, but saveNSx is crashing without a value here
end


for i=1:length(NS5.ElectrodesInfo),
    switch NS5.ElectrodesInfo(i).ConnectorBank
        case 'E'
            NS5.ElectrodesInfo(i).ConnectorBank = 'A';
        case 'F'
            NS5.ElectrodesInfo(i).ConnectorBank = 'B';
        case 'G'
            NS5.ElectrodesInfo(i).ConnectorBank = 'C';
        case 'H'
            NS5.ElectrodesInfo(i).ConnectorBank = 'D';
    end
    NS5.ElectrodesInfo(i).ElectrodeID = NS5.ElectrodesInfo(i).ElectrodeID - 128;
end

outputFileName =  [fullfile(NS5.MetaTags.Filename(1:end)),'-modified',  NS5.MetaTags.FileExt];

outputFilePath = fullfile(fileparts(fileToRead),'_nPlayVisualize/');
if ~exist(outputFilePath,'dir'),
    mkdir(outputFilePath)
end;
outputFilePath = fullfile(outputFilePath,outputFileName);


saveNSx(NS5,outputFilePath);   %saveNSx() currently has a minor bug where the dialog box
%saveNSx(NS5);   %saveNSx() currently has a minor bug where the dialog box
%that asks for a file name doesn't actually do anything.
%The file will save as the original file name with
%'-modified' added before the file extension



end