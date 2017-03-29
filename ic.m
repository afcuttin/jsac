function [outRandomAccessFrame] = ic(raf,sicParameters,receivedBursts)
% function [outRandomAccessFrame,ackedBursts] = sic(raf,sicParameters)
%
% perform Successive Interference Cancellation (SIC) on a given Random Access Frame for Contention Resolution Diversity Slotted Aloha
%
% Input parameters
%   - raf:     the structure of matrices containing slots and bursts informations
%   - maxIter: the maximum number of times the Successive Interference Cancelation can be performed TODO: change the description of the second argument of the function [Issue: https://github.com/afcuttin/crdsa/issues/19]
%
% Output parameters
%   - outRandomAccessFrame: the structure of matrices containing slots and bursts informations, after SIC
%   - ackedBursts.slot:     an array containing the column (slots) indices of acknowledged bursts after SIC
%   - ackedBursts.source:   an array containing the row (sources) indices of acknowledged bursts after SIC

for ii = 1:numel(receivedBursts.slot)
    twinPcktCol = raf.twins{ receivedBursts.source(ii),receivedBursts.slot(ii) };
    for twinPcktIdx = 1:length(twinPcktCol)
        raf.status(receivedBursts.source(ii),twinPcktCol(twinPcktIdx))        = 0; % interference cancelation
        % raf.receivedPower(receivedBursts.source(ii),twinPcktCol(twinPcktIdx)) = sicParameters.residual * raf.receivedPower(receivedBursts.source(ii),twinPcktCol(twinPcktIdx));
        raf.slotStatus(twinPcktCol(twinPcktIdx)) = 2;
    end
end

outRandomAccessFrame = raf;