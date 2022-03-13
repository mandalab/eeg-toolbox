
function [transformedSync] = transformSync(sync,transforms)
% sync = ecogDC12, BRDC12, ecogDC09, BRDC09
% positive value of transforms = trim, negative = NaN-pad
%
%   2/2020 - SJ: changed "<" to ">" in elseif statement: transform_micro(1) > 0
%   6/2020 - SJ: Add kayboard- this function should no longer be called
%

%SJ
fprintf('%s\n','STOP!! If you are here, that means you called transformSync, which should not be called anymore! Please make sure your SVN is up to date and tell SJ!');
keyboard

transform_ecog = transforms{1,1};
transform_micro = transforms{1,2};

if length(find(transform_micro~=0))==0
    transform_micro = -transform_ecog;
    transform_ecog = [0 0 0];
end

micro_dc09 = sync{1,4};
micro_dc12 = sync{1,2};


% transform correction, then drift correction


% NaN-pad or take data out.. straightforward
if transform_micro(1) < 0
    micro_dc09 = cat(2,NaN(1,abs(transform_micro(1))),micro_dc09);
    micro_dc12 = cat(2,NaN(1,abs(transform_micro(1))),micro_dc12);
elseif transform_micro(1) > 0
    micro_dc09 = micro_dc09(transform_micro(1):end);
    micro_dc12 = micro_dc12(transform_micro(1):end);
end
if transform_micro(2) < 0
    micro_dc09 = cat(2,micro_dc09,NaN(1,abs(transform_micro(2))));
    micro_dc12 = cat(2,micro_dc12,NaN(1,abs(transform_micro(2))));
elseif transform_micro(2) > 0
    micro_dc09 = micro_dc09(1:end-transform_micro(2));
    micro_dc12 = micro_dc12(1:end-transform_micro(2));
end

% drift correction.. cue insanity

if transform_micro(3) ~= 0
    driftApplyToMicro = transform_micro(3);
    correctionInterval = floor((length(micro_dc09))/(abs(driftApplyToMicro)));
    %- indices that will get a second copy (steps of correction interval, starting 1/2 way in correction interval
    manipulateHere = zeros(1,(abs(driftApplyToMicro)));
    for i=1:abs(driftApplyToMicro)
        %  put it in the middle of the correction interval because that makes sense
        %       and so an extra sample isn't added to the end on even drifts
        manipulateHere(i) = floor(correctionInterval/2)+1 + correctionInterval*(i-1);
    end
    
    
    % if removing data:
    if driftApplyToMicro < 0
        for sample = 1:length(manipulateHere)
            % lets avg the to-be-removed sample with the next one,
            % and assign that to the next sample
            samplesToChangeAcrossChans = nanmean([micro_dc09(:,manipulateHere(sample)+1),micro_dc09(:,manipulateHere(sample))],2);
            micro_dc09(:,manipulateHere(sample)+1) = samplesToChangeAcrossChans;
            % do same for reref
            samplesToChangeAcrossChans = nanmean([micro_dc12(:,manipulateHere(sample)+1),micro_dc12(:,manipulateHere(sample))],2);
            micro_dc12(:,manipulateHere(sample)+1) = samplesToChangeAcrossChans;
        end
        % now, lets remove all the samples
        micro_dc09(:,manipulateHere) = [];
        micro_dc12(:,manipulateHere) = [];
        
    elseif driftApplyToMicro > 0 %adding data
        
        % prepare new series
        micro_dc09_aligned = NaN(size(micro_dc09,1),size(micro_dc09,2)+abs(driftApplyToMicro));
        micro_dc12_aligned = NaN(size(micro_dc12,1),size(micro_dc12,2)+abs(driftApplyToMicro));
        
        
        for sample = 1:length(manipulateHere)
            % lets avg the to-be-removed sample with the next one,
            % and assign that to the next sample
            samplesToChangeIntoNoreref = nanmean([micro_dc09(:,manipulateHere(sample)+1),micro_dc09(:,manipulateHere(sample))],2);
            % do same for reref
            samplesToChangeIntoProcessed = nanmean([micro_dc12(:,manipulateHere(sample)+1),micro_dc12(:,manipulateHere(sample))],2);

            % now, we need to carefully stitch
            % together real data with new insertions..
            if sample==1
                % insert the new, stretched data into the first slot (larger than the initiial data by 1 samples)
                micro_dc09_aligned(:,1:manipulateHere(sample)+sample) = cat(2,micro_dc09(:,1:manipulateHere(sample)),samplesToChangeIntoNoreref);
                micro_dc12_aligned(:,1:manipulateHere(sample)+sample) = cat(2,micro_dc12(:,1:manipulateHere(sample)),samplesToChangeIntoProcessed);
            elseif sample==length(manipulateHere)
                % if this is last sample, insert end data first
                micro_dc09_aligned(:,manipulateHere(sample)+sample:size(micro_dc09_aligned,2)) = micro_dc09(:,manipulateHere(sample):size(micro_dc09,2));
                micro_dc12_aligned(:,manipulateHere(sample)+sample:size(micro_dc12_aligned,2)) = micro_dc12(:,manipulateHere(sample):size(micro_dc12,2));
                % then insert data up to the last manipulation point, along with newly created sample
                micro_dc09_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,micro_dc09(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoNoreref);
                micro_dc12_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,micro_dc12(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoProcessed);
            else
                % if this is a middle sample, insert it
                % along with newly created sample
                micro_dc09_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,micro_dc09(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoNoreref);
                micro_dc12_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,micro_dc12(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoProcessed);
            end
        end
        
        micro_dc09 = micro_dc09_aligned;
        micro_dc12 = micro_dc12_aligned;
        
    end
    
    
end

transformedSync = {micro_dc09 micro_dc12};

end

