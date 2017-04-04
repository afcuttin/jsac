function [outRandomAccessFrame] = ic(raf,receivedBursts)
% function [outRandomAccessFrame] = ic(raf,receivedBursts)
%
% perform Successive Interference Cancellation (SIC) on a given Random Access Frame
%
% Input parameters
%   - raf:     the structure of matrices containing slots and bursts informations
%   - receivedBursts: list of successfully decoded bursts whose replicas shall be cancelled from the RAF
%
% Output parameters
%   - outRandomAccessFrame: the structure of matrices containing slots and bursts informations, after SIC

for ii = 1:numel(receivedBursts.slot)
    twinPcktCol = raf.twins{ receivedBursts.source(ii),receivedBursts.slot(ii) };
    for twinPcktIdx = 1:length(twinPcktCol)
        raf.status(receivedBursts.source(ii),twinPcktCol(twinPcktIdx))        = 0; % interference cancelation
        raf.slotStatus(twinPcktCol(twinPcktIdx)) = 2;
    end
end

outRandomAccessFrame = raf;