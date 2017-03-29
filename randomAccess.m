function [input] = randomAccess(queuesMatrix)
% function [in] = randomAccess(out)
%
% Simulation of Multiple Random Access
%
% Input parameters:
% * input.sources: number of sources (type: integer)
% * input.queueLength: number of packets (burst) in the fifo queue for each source. The length must be the same for all sources (type: integer)
% * input.maMode: multiple access mode: crdsa or csa (type: string)
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)
% * input.rafLength: the length of the random access frame (RAF) (type: integer)
%
% Output
%
% queuesMatrix: a logical matrix of input.sources rows and input.queueLength columns where each cell is set to 1 if the packet was successfully decoded or 0 if it was not

queueStatus = input.queueLength * ones(input.sources,1);
% queste sono tutte le variabili ereditate dal vecchio codice, se possibile provvedere al refactoring
source.number          = input.sources;
raf.length             = input.rafLength; % Casini et al., 2007, pag.1413
% simulationTime         = 1000; % total number of RAF
% simulationTime         = 10; % for TEST purposes only - comment when doing the real thing
capturePar.status      = 3;
% sicPar.maxIter         = 16;
sicPar.maxIter         = 4;
sicPar.minIter         = 4;
capturePar.threshold   = 14; % the value of the parameter is obtained as follows: thr_val_dB + 11 = capturePar.threshold; thr_val_dB values are -10:1:10
capturePar.criterion   = 'power';
capturePar.type        = 'basic';
% carico il file che contiene le probabilità di cattura
load('capt_SUL.mat','C');
capturePar.probability = C;

for it = sicPar.minIter:sicPar.maxIter

    maxIter = it;

    ackdPacketCount          = 0;
    pcktTransmissionAttempts = 0;
    pcktCollisionCount       = 0;
    source.status            = zeros(1,source.number);
    source.backoff           = zeros(1,source.number);
    % legit source statuses are always non-negative integers and equal to:
    % 0: source has no packet ready to be transmitted (is idle)
    % 1: source has a packet ready to be transmitted, either because new data must be sent or a previously collided packet has waited the backoff time
    % 2: source is backlogged due to previous packets collision
    pcktGenerationTimestamp  = zeros(1,source.number);

    while sum(queueStatus) > 0

    	assert(all(queueStatus >= 0),'The status of a queue cannot be negative');

        raf.status               = zeros(source.number,raf.length); % memoryless
        raf.slotStatus           = int8(zeros(1,raf.length));
        raf.twins                = cell(source.number,raf.length);
        changedSlots             = 0;

        % create the RAF





