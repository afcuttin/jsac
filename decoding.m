function [outRandomAccessFrame,ackedBursts] = decoding(raf,capture)
% function [outputRandomAccessFrame,ackedPcktsCol,ackedPcktsRow] = capture(inputRandomAccessFrame,,capture)
% TODO: update decoding function help [Issue: https://github.com/afcuttin/jsac/issues/14]
%
% perform Successive Interference Cancellation (SIC) on a given Random Access Frame for Contention Resolution Diversity Slotted Aloha
%
% +++ Input parameters
% 		- raf: the structure of matrices containing slots and packets informations
% 		- maxIter: the maximum number of times the Successive Interference Cancelation can be performed
%
% +++ Output parameters
% 		- outRandomAccessFrame: the structure of matrices containing slots and packets informations, after SIC
% 		- ackedPcktsCol: an array containing the column indices of acknowledged packets after SIC
% 		- ackedPcktsRow: an array containing the row indices of acknowledged packets after SIC

ackedBursts.slot   = [];
ackedBursts.source = [];

switch capture.accessMethod
    case 'crdsa'
        for si = 1:raf.length
            % if sum(raf.receivedPower([1:end-1],si)) > 0
            if sum(raf.status(:,si)) >= 1
                collided = find(raf.status(:,si) == 1);
                numCollided = numel(collided);
                assert(numCollided <= size(capture.probability,1),'The number of colliding sources exceeds the maximum value permitted by the scenario, which is %u.',size(capture.probability,1))
                % [~,captured] = max(raf.receivedPower([1:end-1],si));
                captured = collided(randi(numCollided,1));
                assert(numel(captured) == 1,'there should be only one burst captured, there are %u instead',numel(captured));
                % captureRatiodB = 10 * log10(raf.receivedPower(captured,si) / sum(raf.receivedPower([1:end ~= captured],si)));
                % if captureRatiodB >= capture.threshold && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                captureExperiment = rand(1);
                if captureExperiment <= capture.probability(numCollided,capture.threshold) && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                    % update the list of acked bursts
                    ackedBursts.slot   = [ackedBursts.slot,si];
                    ackedBursts.source = [ackedBursts.source,captured];
                    % update the raf
                    raf.status(captured,si)        = 0;
                    % raf.receivedPower(captured,si) = capture.sicResidual * raf.receivedPower(captured,si);
                    % sir has changed, update the slot status
                    raf.slotStatus(si) = 2;
                elseif captureExperiment > capture.probability(numCollided,capture.threshold) || raf.status(captured,si) ~= 1 || ismember(captured,ackedBursts.source)
                    % niente da fare
                    raf.slotStatus(si) = 0;
                else
                    error('Something bad happened');
                end
            else % empty slot, only noise
                % skip this slot
            end
        end
    case 'csa'
        for si = 1:raf.length
            if sum(raf.status(:,si)) >= 1

                % find the number of collided sources
                collided         = find(raf.status(:,si) == 1)
                numCollided      = numel(collided)
                % select one source to be captured and find the slot indices of its corresponding segments
                captured         = collided(randi(numCollided,1))
                assert(numel(captured) == 1,'there should be only one burst captured, there are %u instead',numel(captured));
                twinPcktCol      = raf.twins{ captured,si }
                numberOfSegments = numel(twinPcktCol)+1;
                % evaluate number of colliding sources for every slot where the captured source has a segment
                % raf.status(:,[si twinPcktCol]) % TEST: delete this line after testing
                raf.status(captured,:) % TEST: delete this line after testing
                collided         = [numCollided sum(raf.status(:,twinPcktCol) == 1)];
                assert(numberOfSegments == numel(collided));
                captureExperiment = rand(1);

                if numberOfSegments == 3
                    [~,rateThrInd] = min(abs(capture.rateThrVec - 2));
                    collided % TEST: delete this line after testing
                    rateThrInd % TEST: delete this line after testing
                    if captureExperiment <= capture.probability3seg(collided(1),collided(2),collided(3),rateThrInd) && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                        % update the list of acked bursts
                        ackedBursts.slot   = [ackedBursts.slot,si] % FIXME ripristinare ;
                        ackedBursts.source = [ackedBursts.source,captured] % FIXME ripristinare ;
                        % update the raf
                        raf.status(captured,si)        = 0;
                        % sir has changed, update the slot status
                        raf.slotStatus(si) = 2;
                        raf.slotStatus % TEST: delete this line after testing
                    elseif captureExperiment > capture.probability3seg(collided(1),collided(2),collided(3),rateThrInd) && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                        raf.slotStatus(si) = 0 % FIXME ripristinare ;
                    % elseif captureExperiment > capture.probability(numCollided,capture.threshold) || raf.status(captured,si) ~= 1 || ismember(captured,ackedBursts.source)
                    %     % niente da fare
                    else
                        error('Something bad happened');
                    end
                elseif numberOfSegments == 4

                    [~,rateThrInd] = min(abs(capture.rateThrVec - 2/3));
                    if captureExperiment <= capture.probability4seg(collided(1),collided(2),collided(3),collided(4),rateThrInd) && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                        % update the list of acked bursts
                        ackedBursts.slot   = [ackedBursts.slot,si];
                        ackedBursts.source = [ackedBursts.source,captured];
                        % update the raf
                        raf.status(captured,si)        = 0;
                        % sir has changed, update the slot status
                        raf.slotStatus(si) = 2;
                    elseif captureExperiment > capture.probability4seg(collided(1),collided(2),collided(3),collided(4),rateThrInd) && raf.status(captured,si) == 1 && ~ismember(captured,ackedBursts.source)
                        % cannot decode/capture here any longer
                        raf.slotStatus(si) = 0;
                    % elseif captureExperiment > capture.probability(numCollided,capture.threshold) || raf.status(captured,si) ~= 1 || ismember(captured,ackedBursts.source)
                    %     % niente da fare
                    else
                        error('Something bad happened');
                    end
                else

                    error('There is something wrong with the number of segments');
                end
            else % empty slot, only noise
                % skip this slot
            end
        end
    otherwise
        error('Please select one of the availables access methods (csa, crdsa).');
end

outRandomAccessFrame = raf;
