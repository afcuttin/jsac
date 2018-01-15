function [outRandomAccessFrame,ackedBursts] = decoding(raf,capture)
% function [outRandomAccessFrame,ackedBursts] = capture(raf,capture)
%
% perform Successive Interference Cancellation (SIC) on a given Random Access Frame for Contention Resolution Diversity Slotted Aloha
%
% +++ Input parameters
% 		- raf: the structure of matrices containing slots and packets informations
% 		- capture: the structure containing the parameters of the capture
%
% +++ Output parameters
% 		- outRandomAccessFrame: the structure of matrices containing slots and packets informations, after SIC
% 		- ackedBursts: the structure containing the column (source) and row (slot) indices of acknowledged packets after SIC

ackedBursts.slot   = [];
ackedBursts.source = [];

switch capture.accessMethod
    case 'crdsa'
        % FIXME: raf.slotStatus shall be reset at this point % raf.slotStatus(:) = 0
        for si = 1:raf.length
            if sum(raf.status(:,si)) >= 1
                collided = find(raf.status(:,si) == 1);
                numCollided = numel(collided);
                captured = collided(randi(numCollided,1));
                assert(numel(captured) == 1,'there should be only one burst captured, there are %u instead',numel(captured));
                captureExperiment = rand(1);
                if numCollided <= size(capture.probability,1)
                    assert(numCollided <= size(capture.probability,1),'The number of colliding sources exceeds the maximum value permitted by the scenario, which is %u.',size(capture.probability,1))
                    captureExpThreshold = capture.probability(numCollided,capture.sinrThrInd);
                elseif numCollided > size(capture.probability,1)
                    captureExpThreshold = 0;
                end

                if captureExperiment <= captureExpThreshold && raf.captureAttempts(captured) < 2 && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                    % update the list of acked bursts
                    ackedBursts.slot   = [ackedBursts.slot,si];
                    ackedBursts.source = [ackedBursts.source,captured];
                    % update the raf
                    raf.status(captured,si)        = 0;
                    % sir has changed, update the slot status
                    raf.slotStatus(si) = 2;
                elseif captureExperiment > captureExpThreshold || raf.status(captured,si) ~= 1 || ismember(captured,ackedBursts.source)
                    % niente da fare
                    raf.captureAttempts(captured) = raf.captureAttempts(captured) + 1;
                    raf.slotStatus(si) = 0;
                else
                    % error('Something bad happened');
                end
            else % empty slot, only noise
                % skip this slot
            end
        end
    case 'crdsa-nc' % nc stands for No Capture (original CRDSA access method)
        cleanBurstSlot  = find(sum(raf.status) == 1);

        if numel(cleanBurstSlot) > 0
            raf.slotStatus(cleanBurstSlot) = 1;
            raf.slotStatus(raf.slotStatus == 2) = 0;

            ii = 1;
            while ii <= numel(cleanBurstSlot)
                cleanBurstRow = find(raf.status(:,cleanBurstSlot(ii)));
                % update the list of acked bursts
                ackedBursts.slot   = [ackedBursts.slot,cleanBurstSlot(ii)];
                ackedBursts.source = [ackedBursts.source,cleanBurstRow];
                % update raf
                raf.status(cleanBurstRow,cleanBurstSlot(ii)) = 0;
                % update slot status
                raf.slotStatus(cleanBurstSlot(ii)) = 0;
                % proceed with the possible clean twin removal
                twinPcktCol = raf.twins{ cleanBurstRow,cleanBurstSlot(ii) };
                for twinPcktIdx = 1:length(twinPcktCol)
                    % raf.status(cleanBurstRow,twinPcktCol(twinPcktIdx)) = 0; % interference cancelation
                    if sum(raf.status(:,twinPcktCol(twinPcktIdx))) == 1 % twin burst was a clean burst
                        nonCollTwinInd = find(cleanBurstSlot == twinPcktCol(twinPcktIdx));
                        if ~isempty(nonCollTwinInd)
                            cleanBurstSlot(nonCollTwinInd) = []; %remove the twin burst from the acked bursts list
                        end
                    end
                end
                ii = ii + 1;
            end
            % end
            outRandomAccessFrame = raf;
        elseif numel(cleanBurstSlot) == 0
            % warning('Nothing to do here, exiting')
            raf.slotStatus(:)    = 0;
            outRandomAccessFrame = raf;
        end
    case {'csa-p','csa-pip'} % csa con cattura in parallelo
        for si = 1:size(raf.status,1) % si means "source index" in this case
            if sum(raf.status(si,:)) >= 1
                % find the segments in the slices
                segments         = find(raf.status(si,:) == 1);
                numberOfSegments = numel(segments);
                collided         = zeros(1,numberOfSegments);
                for slin = 1:numberOfSegments
                    collided(slin) = sum(raf.status(:,segments(slin)) == 1);
                end
                % evaluate capture probability
                if numberOfSegments == 3
                    % [~,rateThrInd] = min(abs(capture.rateThrVec - 2));
                    [~,rateThrInd] = min(abs(capture.rateThrVec - 4/3));
                    captureExpThreshold = capture.probability3seg(collided(1),collided(2),collided(3),rateThrInd);
                elseif numberOfSegments == 4
                    % [~,rateThrInd] = min(abs(capture.rateThrVec - 2/3));
                    [~,rateThrInd] = min(abs(capture.rateThrVec - 1));
                    captureExpThreshold = capture.probability4seg(collided(1),collided(2),collided(3),collided(4),rateThrInd);
                else
                    error('There is something wrong with the number of segments, they are %u',numberOfSegments);
                end

                captureExperiment = rand(1);

                if captureExperiment <= captureExpThreshold
                    % update the list of acked bursts
                    ackedBursts.slot   = [ackedBursts.slot,segments(end)];
                    ackedBursts.source = [ackedBursts.source,si];
                    % update the raf
                    raf.status(si,segments) = 0; % NOTE: should be segments(end) as written 3 lines above
                    % sir has changed, update the slot status % NOTE: delete the following two lines, as they are useless in this scenario
                    raf.slotStatus(si) = 2;
                elseif captureExperiment > captureExpThreshold
                    % niente da fare
                    raf.slotStatus(si) = 0; % NOTE: necessario? direi di no
                else
                    error('Something bad happened');
                end
            end
        end
    case 'csa' % csa in serie per la cancellazione iterativa
        raf.slotStatus(:) = 0; % this reset is required to take into account only the cancellations of the current iteration
        % sort the sources in ascending order according to the slice where the last segment of a packet is transmitted
        [srcs,slcs]       = ind2sub(size(raf.status),find(raf.status));
        [srcsSor,srcsInd] = sort(srcs);
        slcsByRow         = slcs(srcsInd); % slices where a segment is present are grouped by rows in the same array
        [srcsDed,srcsOcc] = unique(srcsSor,'last'); % srcsOcc contains the last occurrence (that is, slice) where a segment of the corresponding source is present
        lastSeg           = slcsByRow(srcsOcc);
        [~,lastSegInd]    = sort(lastSeg);
        srcsUns           = srcsSor(srcsOcc);
        srcsSor           = srcsUns(lastSegInd);

        for si = transpose(srcsSor) % si means "source index" in this case
            if sum(raf.status(si,:)) >= 1 % NOTE: this conditional is here to check that the current source has segments in the slices; therefore the condition should check against the (minimum) value of segments that a source can put in the frame
                % find the segments in the slices
                segments         = find(raf.status(si,:) == 1);
                numberOfSegments = numel(segments);
                collided         = zeros(1,numberOfSegments);
                for slin = 1:numberOfSegments
                    collided(slin) = sum(raf.status(:,segments(slin)) == 1);
                end
                % evaluate capture probability
                if numberOfSegments == 3 % NOTE: the value '3' is hard-coded: the code may suffer in case the number of segments is changed
                    collided(collided > size(capture.probability3seg,1)) = size(capture.probability3seg,1); % in case there are more sources than size(capture.probability3seg,1))
                    % [~,rateThrInd] = min(abs(capture.rateThrVec - 2)); % NOTE: superseded by csa with qpsk
                    % captureExpThreshold = capture.probability3seg(collided(1),collided(2),collided(3),rateThrInd); % NOTE: superseded by csa with qpsk
                    value1=capture.probability3seg(collided(1),collided(2),collided(3),capture.rt_ind(3));
                    value2=capture.probability3seg(collided(1),collided(2),collided(3),capture.rt_ind(3)+1);
                    captureExpThreshold = value1+(raf.numSegments*capture.bit_per_symbol/3-capture.rateThrVec(capture.rt_ind(3)))*(value2-value1)/(capture.rateThrVec(capture.rt_ind(3)+1)-capture.rateThrVec(capture.rt_ind(3)));
                elseif numberOfSegments == 4 % NOTE: the value '4' is hard-coded: the code may suffer in case the number of segments is changed
                    collided(collided > size(capture.probability4seg,1)) = size(capture.probability4seg,1); % in case there are more sources than size(capture.probability4seg,1))
                    % [~,rateThrInd] = min(abs(capture.rateThrVec - 2/3)); % NOTE: superseded by csa with qpsk
                    % captureExpThreshold = capture.probability4seg(collided(1),collided(2),collided(3),collided(4),rateThrInd); % NOTE: superseded by csa with qpsk
                    value1=capture.probability4seg(collided(1),collided(2),collided(3),collided(4),capture.rt_ind(4));
                    value2=capture.probability4seg(collided(1),collided(2),collided(3),collided(4),capture.rt_ind(4)+1);
                    captureExpThreshold=value1+(raf.numSegments*capture.bit_per_symbol/4-capture.rateThrVec(capture.rt_ind(4)))*(value2-value1)/(capture.rateThrVec(capture.rt_ind(4)+1)-capture.rateThrVec(capture.rt_ind(4)));
                else
                    error('There is something wrong with the number of segments, they are %u',numberOfSegments);
                end

                captureExperiment = rand(1);

                if captureExperiment <= captureExpThreshold && raf.captureAttempts(si) < 2
                    % update the list of acked bursts
                    ackedBursts.slot   = [ackedBursts.slot,segments(end)];
                    ackedBursts.source = [ackedBursts.source,si];
                    % update the raf
                    raf.status(si,segments) = 0; % NOTE: should be segments(end) as written 3 lines above
                    % sir has changed, update the slot status % NOTE: delete the following two lines, as they are useless in this scenario
                    raf.slotStatus(si) = 2;
                elseif captureExperiment > captureExpThreshold
                    % niente da fare
                    raf.captureAttempts(si) = raf.captureAttempts(si) + 1;
                    raf.slotStatus(si) = 0; % NOTE: necessario? direi di no
                else
                    error('Something bad happened');
                end
            end
        end
    case 'csa-nc' % nc stands for No Capture (original CSA access method)
        [srcs,slcs]        = ind2sub(size(raf.status),find(raf.status));
        [srcsSor,srcsInd]  = sort(srcs);
        slcsByRow          = slcs(srcsInd); % slices where a segment is present are grouped by rows in the same array
        clnSlcs            = find(sum(raf.status) == 1); % find clean slices
        clnSrcsRaw         = srcsSor .* ismember(slcsByRow,clnSlcs); % find sources associated with clean slices
        clnSrcs            = clnSrcsRaw(clnSrcsRaw ~= 0);
        [counts,countsIdx] = histc(clnSrcs, unique(clnSrcs)); % evaluate number of clean slices for each source
        ackedBursts.source = transpose(unique(clnSrcs(counts(countsIdx) >= raf.numSegments))); % successfull sources in column
        [~,~,slcsIdx]      = intersect(ackedBursts.source,clnSrcsRaw);
        ackedBursts.slot   = transpose(slcsByRow(slcsIdx));
        % update raf
        if ~isempty(ackedBursts.source)
            raf.status(sub2ind(size(raf.status),ackedBursts.source,ackedBursts.slot)) = 0;
        end
        % reset slot statuses
        raf.slotStatus(:) = 0;
        outRandomAccessFrame = raf;
    otherwise
        error('Please select one of the availables access methods (csa, crdsa).');
end
outRandomAccessFrame = raf;
