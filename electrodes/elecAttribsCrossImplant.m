function [varargout] = elecAttribsCrossImplant(subj, rootEEGdir, srcImplant, D, attributes, silent)
% elecAttribsCrossImplant is designed for transferring
% attributes (isResected/ictal/interictal) from one implant to another.
% For example, if the attributes are known only for the 2nd implant and they need
% to be known for the first implant.
%
% Attributes for the unknown implant are inferred from the known implant as follows:
%   If electrode E is known to have attribute A at location XYZ during one implant,
%   then all electrodes within DIST of XYZ are assigned attribute A for the other,
%   unknown umplant.
%
% USAGE:
%   elecAttribsCrossImplant(subj, rootEEGdir, srcAsOf, [dist=6])
%
% INPUT:
%   subj
%   rootEEGdir
%   srcImplant - 1 or 2. If 1, the first implant is the source of known attributes and
%       we want to estimate the 2nd implant
%   D - Optional (default 6). Distance in mm to use for new attributes
%
%   attributes - Optional. Default {'resected','interictal','ictal'}. Limit results to these
%               attributes
%   silent    - If 1, do not print any output. Default false
%
% OUTPUT
%   srcChan     - n length cell array of chan names for source of attribute
%   cp2Chan     - n length cell array of chan names to be copied with attribute
%   attribVec   - n length cell array of 3-length logical vectors giving the attributes as:
%                 [ictal, interictal, resected]. e.g. [1 0 1] means source was ictal and resected
% REVISION HISTORY:
%   05/18 MST - Created

DEBUG = 0;
DATE_FORMAT = 'MM/dd/yy';

if nargin < 4 || isempty(D)
    D = 6;
end
if nargin < 5 || isempty(attributes)
    attributes = {'resected','interictal','ictal'};
elseif ischar(attributes)
    attributes = cellstr(attributes);
end
assert(numel(attributes) > 0, 'Must pass cell array subset of resected, interictal, ictal for attributes');
if nargin < 6 || isempty(silent)
    silent = 0;
end

assert(srcImplant == 1 || srcImplant == 2, 'srcImplant must be 1 or 2');
info = getElementInfo(subj, rootEEGdir);
assert(~isempty(info), 'Empty element Info');




%% Divide hardware into source and copyTo sets
all_dates = union(info.dateIn, info.dateOut);
all_dts = datetime(all_dates, 'InputFormat', DATE_FORMAT);
all_dts = sort(all_dts(~isnat(all_dts)));
assert(numel(all_dts) >= 2, 'There must be at least 2 dates in elementInfo for there to have been 2 implants');
assert(numel(all_dts) <= 3, 'There are more than three dates in elementInfo; this should not be the case for 2 implants');

dateIn_dt  = datetime(info.dateIn,  'InputFormat', DATE_FORMAT);
dateOut_dt = datetime(info.dateOut, 'InputFormat', DATE_FORMAT);

phys_mask = strcmpi(info.chanType, 'PHYS');

% first implant defined as "in on first date, out on 2nd date"
im1_mask = (dateIn_dt == all_dts(1)) & (dateOut_dt == all_dts(2));

% 2nd implant defined as "in after the first date"
im2_mask = (dateIn_dt > all_dts(1));

% present in Both implants defined as "in on first date and out after the 2nd date"
imBoth_mask = (dateIn_dt == all_dts(1)) & (isnat(dateOut_dt) | (dateOut_dt > all_dts(2)));

sh_mask = contains(info.tagName, '_sh');
sh_tags = cellfun(@(x) strrep(x,'_sh',''), info.tagName(sh_mask), 'uniformOutput',0);
pre_sh_mask = ismember(info.tagName, sh_tags);

if DEBUG
    print(silent, 'Implant 1:\n');
    info(phys_mask & im1_mask, :)
    print(silent, 'Implant both:\n');
    info(phys_mask & imBoth_mask,:)
    print(silent, 'Implant 2:\n');
    info(phys_mask & im2_mask,:)
    disp('sh tags');
    disp(info.tagName(sh_mask));
    disp('Pre-shift tags');
    disp(info.tagName(pre_sh_mask));
end


if srcImplant == 1
    src_mask = (im1_mask | imBoth_mask) & phys_mask;
    cp2_mask = im2_mask & phys_mask;
    cp2Implant = 2;
else
    src_mask = (im2_mask | imBoth_mask) & phys_mask;
    cp2Implant = 1;
    % If implant-2 is source, we should restrict copyTo electrodes (implant-1) to channels
    % that were shifted. For example, if channel X was in implant 1 but was never
    % shifted, it would not have been missed during manual attribute assignment.
    % Whereas if it was shifted, the manual annotator should mark X_sh and thus
    % we really do need to infer X's attribute.
    cp2_mask = im1_mask & phys_mask & pre_sh_mask;
end

src_tag = info.tagName(src_mask);
cp2_tag = info.tagName(cp2_mask);


if DEBUG
    disp('source tags');
    disp(src_tag);
    disp('copyTo tags');
    disp(cp2_tag);
end


%% Get attributes
jack = getJackTable(subj, rootEEGdir);
leadsIctal = getLeads(subj, rootEEGdir, 'jackTable',jack, 'markedAs','isIctal');
leadsInter = getLeads(subj, rootEEGdir, 'jackTable',jack, 'markedAs','isInterictal');
leadsResec = getLeads(subj, rootEEGdir, 'jackTable',jack, 'markedAs','isResected');

%% Get source XYZ
bd = braindata2(subj, rootEEGdir, 1);
if isempty(bd.tal), 
    fprintf('\n error, cant push resection/ictal/interictal data from last implant to prior implant without tal being processed');
    keyboard;
    return; 
end
t_src = bd.tal.xyz(ismember(util_split_stringnum(bd.tal.xyz.chanName), src_tag), :);
t_cp2 = bd.tal.xyz(ismember(util_split_stringnum(bd.tal.xyz.chanName), cp2_tag), :);

% Normally, if an electrode is present in both implants, the 2nd implant's XYZ takes precedence
% For this case, if implant 1 is source, we want its xyz to take precedence. This is stored in braindata
if srcImplant == 1
    if isprop(bd, 'zloc') && isfield(bd.zloc, 'proj_1_xyz')
        t1_src = bd.zloc.proj_1_xyz(ismember(util_split_stringnum(bd.zloc.proj_1_xyz.chanName), src_tag), :);
        for i = 1:height(t1_src)
            replace_row_ndx = find(strcmpi(t_src.chanName, t1_src.chanName{i}), 1);
            t_src(replace_row_ndx, {'chanName','x','y','z'}) = t1_src(i, {'chanName','x','y','z'});
        end
    else
        warning('Looking for proj_1_xyz in braindata and not found! Possibly using overwritten xyz locations');
    end
end


xyz_src = t_src{:, {'x','y','z'}};
xyz_cp2 = t_cp2{:, {'x','y','z'}};
cross_dists = pdist2(xyz_src, xyz_cp2);

%% Apply attributes to all copyTo electrodes within D of each source electrode
[ndxs_src,ndxs_cp2] = find(cross_dists <= D);
if DEBUG
    print(silent, '%d total electrodes found within %d mm\n', numel(ndxs_src), D);
end

print(silent, '-----------------------------------------------\n');
print(silent, 'The following electrodes from implant %d have attributes\n', srcImplant);
print(silent, 'that should be transfered to implant %d and are within %d mm:\n', cp2Implant, D);
print(silent, '-----------------------------------------------\n');
print(silent, 'SRC --> COPYTO (DIST) ICTAL INTERICTAL RESECTED\n\n');
cnt = 0;
attrib_requested_vec = [ismember('ictal', attributes), ismember('interictal', attributes), ismember('resected', attributes)];
srcChan = {}; %- JW converted to empty cell arays (was initalized as empty vectors, which caused problems downstream if expecting to compare cells)
cp2Chan = {};
attribVec = {};
for k = 1 : numel(ndxs_src)
    i_src = ndxs_src(k);
    j_cp2 = ndxs_cp2(k);
    chan_src = t_src.chanName{i_src};
    chan_cp2 = t_cp2.chanName{j_cp2};
    isIctal = ismember(chan_src, leadsIctal);
    isInter = ismember(chan_src, leadsInter);
    isResec = ismember(chan_src, leadsResec);
    attrib_present_vec = [isIctal, isInter, isResec];
    if any(dot(double(attrib_present_vec), double(attrib_requested_vec)))
        cnt = cnt + 1;
        
        %- JW addes this text output when just resections are being requested (as is done in eegPrepAndAlign) so easier for user to figure out what is missing
        strDestIsResect = '';
        if length(attributes)==1 & strcmp(attributes{1},'resected') & ismember(chan_cp2, leadsResec),
            strDestIsResect = '  << (resection status already copied to destination)';
        end
        
        print(silent, '%8s --> %8s (%0.1f mm) %d %d %d%s\n', ...
            chan_src, chan_cp2, cross_dists(i_src,j_cp2), isIctal, isInter, isResec,strDestIsResect);
        
        srcChan{cnt} = chan_src;
        cp2Chan{cnt} = chan_cp2;
        attribVec{cnt} = attrib_present_vec;
    end
end
if cnt == 0
    print(silent, 'No electrodes found!\n');
end
print(silent, '-----------------------------------------------\n');

if nargout >= 1, varargout{1} = srcChan; end
if nargout >= 2, varargout{2} = cp2Chan; end
if nargout >= 3, varargout{3} = attribVec; end

end

function print(silent, varargin);
if ~silent
    fprintf(varargin{:});
end
end

function other_name = getShiftname(subj, rootEEGdir, chan, implant)
% get the name of chan in implant
% example: chan=G, implant=2, other_name may be G_sh
%          chan=G, implant=1, other_name will be G
filename = fullfile(rootEEGdir, subj, 'tal/zloc/mr_pre/shift_rename.csv');
if ~exist(filename, 'file')
    other_name = chan;
    return
end
t = readtable(filename);
irow = find(strcmpi(t(:,1), chan), 1);
if isempty(irow)
    irow = find(strcmpi(t(:,2), chan), 1);
end
if isempty(irow)
    other_name = chan;
    return
end
if implant == 1
    other_name = t{irow, 2};
else
    other_name = t{irow, 1};
end
end
