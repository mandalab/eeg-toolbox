function auto_annotate_old(wavdir)
    % Author: Robert Yafe
    % comments and revisions: Mike
    
    fprintf('Running auto_annotate_old: **This is for NIH030 and NIH032 ONLY**\n');
    % In this version, we want to create an .ann for every .wav
    % Every .wav contains a single trial
    
    % in this version, a single trial looks like this:
%     1437747816361	0	FIXATION_ON
%     1437747816878	0	FIXATION_OFF
%     1437747817394	0	PROBEWORD_ON	CLOCK	TARGET_MATCH	BRICK	0_000_10(01)_14377478163
%     1437747817398	1	REC_START
%     1437747819061	1	REC_END
%     1437747819079	0	MATCHWORD_ON	BRICK
%     1437747819395	0	PROBEWORD_OFF	CLOCK	MATCHWORD_OFF	BRICK
%     1437747819911	0	FIXATION_ON

% so in our line loop, we want to store TARGET_MATCH and ann file (words end-1 and end) on the probeword_on line
    % record the rec_start
    % record the rec_end
    % then annotate

    
    % ** define minh here if you want script to run without asking **
    minh = 2;
    MARK_ALL_INCOMPLETE = 0;
    RERUN_ALL = 1; % delete all .ann and .tmp
    DBG = 0; % if 1, show figures and print session.log 
    % ** ^ comment this out to have script ask you for it **

    words = {'BRICK','CLOCK','GLASS','JUICE','PANTS'};
    
    if ~exist('wavdir','var')
        %wavdir = '/Volumes/shares/FRNU/dataWorking/eeg/needsAnnotation/paRemap annotation/NIH037/session_4/';
        %wavdir = '/Volumes/Shares/FRNU/dataWorking/eeg/old/NIH043/behavioral/paRemap/session_1/';
        %wavdir = '~/NIH/dataWorking/eeg/NIH037/behavioral/paRemap/session_0/';
        
    end

    % findpeaks is overloaded in eeg_toolbox, so grab it directly from signal toolbox
    cur = pwd;
    cd(fullfile(toolboxdir('signal'), 'signal'));
    findpeaks_signal = @findpeaks;
    cd(cur);
    
    
    fprintf('********************\n')
    fprintf('wavdir is: %s\n', wavdir);
    fprintf('Exists: %d\n', exist(wavdir,'dir') > 0);
    if exist('minh', 'var'), fprintf('minh set to: %d\n', minh); end
    if ~exist(wavdir,'dir'), return; end
    fprintf('********************\n')
    
    if wavdir(end) ~= '/'
        wavdir = [ wavdir '/' ];
    end
    
    logfile = fullfile(wavdir, 'session.log');
    f = fopen(logfile);
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

    
    % For every line
    while ischar(l)
        
        if DBG
            disp(l);
        end
        
        % Assume format is [timestamp, integer, 2 strings]
        xTOT = textscan(l,'%f%d%s%s');
        type = xTOT{3}{1};
        mstime = xTOT{1}(1);
        
        
        
%         % store word list
%         if strcmp(type, 'WORD_LIST')
%             w = strsplit(l);
%             words = w(end-3:end);
%         end
        
        % start wav recording
        if strcmp(type,'REC_START')
            
            % store recording window start
            made_annotation = 0;
            
            % last line was important, grab the annotation file and
            % the target word now. Also probe_time
            % (last line is % timestamp	0	PROBEWORD_ON	CLOCK	TARGET_MATCH	BRICK	filename)
            if isempty(xTOT{4})
                wlast = strsplit(llast);
                xTOT{4} = wlast(end);
                target_word = wlast{end-1};
                probe_time = str2double(wlast{1});
            end
            
            recfile = [wavdir xTOT{4}{1} '.wav'];
            tmpfile = [wavdir xTOT{4}{1} '.tmp'];
            annfile = [wavdir xTOT{4}{1} '.ann'];
            write_to_file = 0;
            
            if RERUN_ALL
                if exist(tmpfile,'file'), delete(tmpfile); end
                if exist(annfile,'file'), delete(annfile); end
            end
            
            % never annotated this before: open an annotation and 
            % print a default header
            if ~exist(tmpfile,'file') && ~exist(annfile,'file')
                f_ann = fopen(annfile,'w');
                for xx = 1:8
                    fprintf(f_ann,'%s\n',hdr{xx});
                end
                fprintf(f_ann,'%s\n','');
                write_to_file = 1;
            end
            
            
            if MARK_ALL_INCOMPLETE
                mark_incomplete(annfile);
            end
            
            
            % read the audio and make markings
            [y,fs] = audioread(recfile);
            t = (1:numel(y))/fs;
            t = t';
            infile = 1;
            z = zscore(y);
            if DBG
                figure(1)
                clf
                plot(t,z);
                hold on
            end
            x = windowed_rms(z,numperwin,numperslide);
            
            t2=numperwin:numperslide:length(y);
            %t2 = (t2-numperwin/2)/fs;
            t2 = t2/fs;

            if DBG
                figure(2)
                clf
                plot(t2,x)
                hold on
                line(get(gca,'xlim'),[2 2],'Color','m')
                line(get(gca,'xlim'),[2 2],'Color','r')
                line(get(gca,'xlim'),[3 3],'Color','g')
                line(get(gca,'xlim'),[4 4],'Color','k')
                pause(.1)
            end

            if ~exist('minh','var')
                minh = input('Enter height threshold: ');
            end
            
            [~,locs] = findpeaks_signal(x,'MinPeakHeight',minh);
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
            
            if isempty(markings)
                % redo with plot and prompt
                figure(2)
                clf
                plot(t2,x)
                hold on
                line(get(gca,'xlim'),[2 2],'Color','m')
                line(get(gca,'xlim'),[2 2],'Color','r')
                line(get(gca,'xlim'),[3 3],'Color','g')
                line(get(gca,'xlim'),[4 4],'Color','k')
                pause(.1)
                
                minh = input('Enter height threshold: ');
                
                [~,locs] = findpeaks_signal(x,'MinPeakHeight',minh);
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
            
            
            
        end % end rec_start line 

        % enter this clause after the inital rec_start and on all subsequent? lines
        % until forever?
        % ? because maybe Mike should change this logic..
        
        if infile
            
            % ms_start is the recording time relative to cue onset (probe_time)
            if strcmp(type,'REC_START')
                rec_offset_ms = mstime - probe_time;
                markings = markings + rec_offset_ms; % make it relative to cue
            end


            % enter this loop once, only on each rec_end line
            if strcmp(type,'REC_END')
                has_bounds = 0;
                % in this old case, the whole audio occurs during the correct time, so find all markings
                mark_inds = find(markings);
                
                % if a mark occurs
                if ~isempty(mark_inds)
                    mark = markings(mark_inds(1));
                    if write_to_file
                        made_annotation = 1;
                        fprintf(f_ann,'%s\t%d\t%s\n',num2str(mark),find(ismember(words,target_word)),target_word);
                        mark_incomplete(annfile); % the user needs to double check this
                    end
                end
                
                % Close annotation .ann file
                if strcmp(type,'REC_END')
                    if ~made_annotation
                        keyboard;
                    end
                    if exist('f_ann','var')
                        try fclose(f_ann); catch, end
                    end
                end
                
            end
        end
        llast = l;
        l = fgetl(f);
    end
    fclose(f);
    
    fprintf('Complete!\n');
end

function mark_incomplete(filename)
    % After an initial edit is made to an annotation, changes are stored in a
    % tmp file until 'Mark Complete'. So to make files appear not complete,
    % create a temp file by duplicating the .ann file
    
    [d,name,ext] = fileparts(filename);
    filename_temp = [fullfile(d, name) '.tmp'];
    if exist(filename, 'file') && ~exist(filename_temp, 'file')
        movefile(filename, filename_temp);
    end

end


