function auto_annotate(wavdir)
% auto_annotate(wavdir)
    
    % ** define minh here if you want script to run without asking **
    %minh = 2.00;
    MARK_ALL_INCOMPLETE = 0;
    RERUN_ALL = 0; % delete all .ann and .tmp
    H_PER_FILE = 0;
    DBG = 0;
    % ** ^ comment this out to have script ask you for it **

    if ~exist('wavdir','var')
        %wavdir = '/Volumes/shares/FRNU/dataWorking/eeg/needsAnnotation/paRemap annotation/NIH037/session_4/';
        wavdir = '/Volumes/Shares/FRNU/dataWorking/eeg/NIH042/behavioral/paRemap/session_0/';
        %wavdir = '~/NIH/dataWorking/eeg/NIH037/behavioral/paRemap/session_0/';
    end
    
    words = {'BRICK','CLOCK','GLASS','JUICE','PANTS'};

%     % findpeaks is overloaded in eeg_toolbox, so grab it directly from signal toolbox
%     cur = pwd;
%     cd(fullfile(toolboxdir('signal'), 'signal'));
%     findpeaks_signal = @findpeaks;
%     cd(cur);
    % *** this issue was resolved by renaming the eeg_toolbox version "findpeaksEEG"
    
    fprintf('********************\n')
    fprintf('wavdir is: %s\n', wavdir);
    fprintf('Exists: %d\n', exist(wavdir,'dir') > 0);
    if exist('minh', 'var'), fprintf('minh set to: %d\n', minh); end
    if ~exist(wavdir,'dir'), return; end
    fprintf('********************\n')
    
    logfile = fullfile(wavdir, 'session.log');
    f = fopen(logfile,'r');
    l = fgetl(f);
    infile = 0;
    has_bounds = 0;

    numperwin = 1000;
    numperslide = 100;


    hdr{1} = '#Begin Header. [Do not edit before this line. Never edit with an instance of the program open.]';
    hdr{2} = '#Annotator:  robot';
    hdr{3} = '#UTC Locally Formatted: July 19, 2016	2:01:53 PM EDT';
    hdr{4} = '#UNIX: 1468951313';
    hdr{5} = '#Program Version: 1.10';
    hdr{6} = '#OS Properties: {os.name	Mac OS X}{os.arch	x86_64}{user.name	auto-annotate}{user.country	US}{user.language	en}';
    hdr{7} = '#Java Properties: {java.runtime.version	1.6.0_65-b14-466.1-11M4716}{java.specification.name	Java Platform API Specification}{java.specification.vendor	Sun Microsystems Inc.}{java.specification.version	1.6}{java.vendor	Apple Inc.}{java.version	1.6.0_65}{java.vm.name	Java HotSpot(TM) 64-Bit Server VM}{java.vm.specification.name	Java Virtual Machine Specification}{java.vm.specification.vendor	Sun Microsystems Inc.}{java.vm.specification.version	1.0}{java.vm.vendor	Apple Inc.}{java.vm.version	20.65-b04-466.1}';
    hdr{8} = '#83 112 97 110 58 32 49 52 54 56 57 53 49 50 57 56 57 50 53 45 49 52 54 56 57 53 49 53 54 56 54 52 56';

    while ischar(l)
        xTOT = textscan(l,'%f%d%s%s');
        type = xTOT{3}{1};
        mstime = xTOT{1}(1);
        if strcmp(type,'REC_END')
            if exist('f_ann','var')
                try fclose(f_ann); catch, end
            end
        end
        if strcmp(type,'REC_START')
            recstarttime = mstime;
            recfile = fullfile(wavdir, [xTOT{4}{1} '.wav']);
            tmpfile = fullfile(wavdir, [xTOT{4}{1} '.tmp']);
            annfile = fullfile(wavdir, [xTOT{4}{1} '.ann']);
            write_to_file = 0;
            if ~exist(tmpfile,'file') && ~exist(annfile,'file')
                f_ann = fopen(annfile,'w');
                for xx = 1:8
                    fprintf(f_ann,'%s\n',hdr{xx});
                end
                fprintf(f_ann,'%s\n','');
                write_to_file = 1;
            elseif RERUN_ALL
                if exist(tmpfile,'file'), delete(tmpfile); end
                if exist(annfile,'file'), delete(annfile); end
            else
                l = fgetl(f);
                continue;
            end
            
            
            if MARK_ALL_INCOMPLETE
                mark_incomplete(annfile);
            end
            
            
            
            [y,fs] = audioread(recfile);
            t = (1:numel(y))/fs;
            t = t';
            t2=numperwin:numperslide:length(y);
            %t2 = (t2-numperwin/2)/fs;
            t2 = t2/fs;
            
            infile = 1;
            z = zscore(y);
            x = windowed_rms(z,numperwin,numperslide);
            
            % remove saturation
%             yy = y;
%             Z_SAT = 10;
%             while max(z) >= Z_SAT
%                 % replace all elements over 10 sigma with random elements
%                 ind_sat = find(z >= Z_SAT);
%                 ysubsample = y((z < Z_SAT) & (z > 0));
%                 sample = ysubsample(randi(length(ysubsample), length(ind_sat), 1));
%                 yy(ind_sat) = sample;
%                 z = zscore(yy);
%                 
%                 figure(2)
%                 clf
%                 plot(t2,x)
%                 hold on
%                 line(get(gca,'xlim'),[2 2],'Color','m')
%                 line(get(gca,'xlim'),[2 2],'Color','r')
%                 line(get(gca,'xlim'),[3 3],'Color','g')
%                 line(get(gca,'xlim'),[4 4],'Color','k')
%                 pause(0.1);
%             end
            
            if DBG
                figure(1)
                clf
                plot(t,z);
                hold on
            end
            

            if DBG
                figure(2)
                clf
                plot(t2,x)
                hold on
                line(get(gca,'xlim'),[2 2],'Color','m')
                line(get(gca,'xlim'),[2 2],'Color','r')
                line(get(gca,'xlim'),[3 3],'Color','g')
                line(get(gca,'xlim'),[4 4],'Color','k')
                set(gca,'YLim',[0,6]);
                pause(.1)
            end

            if ~exist('minh','var') || H_PER_FILE
                minh = input('Enter height threshold: ');
                
            end
            
            [~,locs] = findpeaks(x,'MinPeakHeight',minh);
            th = prctile(x,60);

            markings = nan(numel(locs),1);
            for i = 1:numel(locs)
                loc = locs(i);
                q = 1;
                while q>th
                    loc = loc-1;
                    if loc == 0
                        q=th;
                        loc = 1;
                    else
                        q = x(loc);
                    end
                end
                markings(i) = t2(loc);
            end
            markings = unique(markings);
            markings = markings*1000;
        end

        if infile
            if strcmp(type,'PROBEWORD_ON')
                ms_start = mstime-recstarttime;
            elseif strcmp(type,'MATCHWORD_ON')
                word = xTOT{4}{1};
                ms_end = mstime-recstarttime;
                has_bounds = 1;
            end
            if has_bounds
                has_bounds = 0;
                mark_inds = find(markings>ms_start & markings<ms_end);
                if ~isempty(mark_inds)
                    mark = markings(mark_inds(1));
                    if write_to_file
                        fprintf(f_ann,'%s\t%d\t%s\n',num2str(mark),find(ismember(words,word)),word);
                        mark_incomplete(annfile);
                    end
                end
            end
        end
        l = fgetl(f);
    end
    fclose(f);
end

function mark_incomplete(filename)
    % After an initial edit is made to an annotation, changes are stored in a
    % tmp file until 'Mark Complete'. So to make files appear not complete,
    % create a temp file by duplicating the .ann file
    
    [d,name] = fileparts(filename);
    filename_temp = [fullfile(d, name) '.tmp'];
    if exist(filename, 'file') && ~exist(filename_temp, 'file')
        movefile(filename, filename_temp);
    end

end

function result = windowed_rms(data,points_per_win,points_per_slide)
    %This computes the moving window average of data which is time X
    %channels
    total_points = size(data,1);
    n_channels = size(data,2);
    points_per_overlap = points_per_win - points_per_slide;
    n_windows = floor((total_points - points_per_overlap)/points_per_slide);
    result = zeros(n_windows,n_channels);
    for i = 1:n_windows
        start = (i-1)*points_per_slide + 1;
        result(i,:) = rms_norm(data(start:start+points_per_win-1,:));
    end
end

function [r,x] = rms_norm(x)
    %takes the root mean square of a vector or a matrix along the first
    %dimension
    meanx = mean(x,1);
    x = x - repmat(meanx,size(x,1),1);
    r = sqrt(mean(x.^2,1));
end