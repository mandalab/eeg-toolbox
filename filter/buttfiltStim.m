function [dat,filters] = buttfiltStim(dat,freqrange,samplerate,filttype,order,datZeroLogic)
%  jwNotch Stim
%
%
%  apply a butterworth filter to stimulation data

%   intended to filter out line noise with stimulation data
%   acts like "buttfilt", but zeros out the stim artifact before filtering
%    that way the ringing caused by the jump to zero is contained within the zero
%
%   Note: tried a version where filts in one direction, zeros again, flips, filts, then flips, but no improvement for high Freqs
%
%   Input: single time series (dat = Ntime x 1 vector) or matrix of time series (dat = Ntime x Nchan)
%
%
%
%BUTTFILT - Wrapper to Butterworth filter.
%
% Butterworth filter wrapper function with zero phase distortion.
%
% FUNCTION:
%   y = buttfilt(dat,freqrange,samplerate,filttype,order)
%   y = buttfilt(dat,filters)
%
% INPUT ARGS: (defaults shown):
%   dat = dat;                % data to be filtered (if data is a matrix, BUTTFILT filters across rows)
%   freqrange = [58 62];      % filter range (depends on type)
%   samplerate = 256;         % sampling frequency
%   filttype = 'stop';        % type of filter ('bandpass','low','high','stop') 
%   order = 4;                % order of the butterworth filter
%   filters = {B,A};          % a set of filters to use on data (created by a previous buttfilt call)
%
%   datZeroLogic = empty if no zeroing time point, else a logical index into dat to zero
%
% OUTPUT ARGS::
%   y = the filtered data
%

% 12/1/04 - PBS - Will now filter multiple times if requested

if ~exist('filttype','var') || isempty(filttype)
    filttype = 'stop';
elseif strcmp(filttype, 'bandpass') && isvector(freqrange)
    if freqrange(1) <= 0
        filttype = 'low';
        freqrange = freqrange(2);
    elseif isinf(freqrange(2))
        filttype = 'high';
        freqrange = freqrange(1);
    end
end

if ~exist('order','var') || isempty(order)
    order = 4;
end

if ~exist('datZeroLogic','var') || isempty(datZeroLogic),
    datZeroLogic = [];
end

if( ~iscell(freqrange) )
    % premade filters were not passed in, so make filters
    nyq = samplerate/2; % Nyquist frequency

    filters = cell(size(freqrange,1),1);
    for i = 1:size(freqrange,1)
        [filters{i,1}, filters{i,2}] = butter(order, freqrange(i,:)/nyq, filttype);
    end
else
    % premade filters were passed in as second argument, so use them
    filters = freqrange;
end


%- data check
if sum(isnan(dat(:)))>0, 
    fprintf('\n data has nan... cant filter with this code'); keyboard; return; 
end


% run the filtfilt for zero phase distortion
for i = 1:size(filters,1)
    
    % following code works for single time series (dat = Ntime x 1 vector) or matrix of time series (dat = Ntime x Nchan)
    if isempty(datZeroLogic) | sum(datZeroLogic)==0,
        dat = filtfilt(filters{i,1},filters{i,2},dat);
    else
        
        %- just zero out stim, filter, then put back what was cut.  Attempt below with separate forwared and reverse didn't make it any better.  Reduce ringing duration at steps with a winder filter range (50-70 instead of 58-62)
        zeroDat = dat(datZeroLogic==1,:);
        dat(datZeroLogic==1,:)=0;
        dat = filtfilt(filters{i,1},filters{i,2},dat);
        dat(datZeroLogic==1,:)=zeroDat;
        
        %- zero stim; filter in one direction; zero stim; flip; filter; flip; return.... didn't improve anything
        % dat(datZeroLogic==1) = 0;
        %dat = filter(filters{i,1},filters{i,2},dat(end:-1:1));
        %dat = dat(end:-1:1);
        %dat(datZeroLogic==1) = 0; %%  was this done incorrectly??... should have flipped logic too?; no, looks ok. dat was reflipped first
        % dat = filter(filters{i,1},filters{i,2},dat);
        
    end
end

