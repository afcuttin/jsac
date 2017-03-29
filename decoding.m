function [outRandomAccessFrame,ackedBursts] = decoding(raf,capture)
% function [outputRandomAccessFrame,ackedPcktsCol,ackedPcktsRow] = capture(inputRandomAccessFrame,,capture)
% TODO: update decoding function help
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

for si = 1:raf.length
    % if sum(raf.receivedPower([1:end-1],si)) > 0
    if sum(raf.status(:,si)) >= 1
        collided = find(raf.status(:,si) == 1);
        numCollided = numel(collided);
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
            % raf.receivedPower(captured,si) = capture.sicResidual * raf.receivedPower(captured,si); % TODO: il pacchetto catturato che residuo lascia? leggere zanella-zorzi
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

outRandomAccessFrame = raf;