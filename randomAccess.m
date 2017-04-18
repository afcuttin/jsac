function [output] = randomAccess(numberOfSources,queueLength,linkMode)
% function [in] = randomAccess(out)
%
% Simulation of Multiple Random Access
%
% Input parameters:
% * numberOfSources: number of sources (type: integer)
% * queueLength: number of packets (burst) in the fifo queue for each source. The length must be the same for all sources (type: integer). Can be a scalar (all the queues have the same length) or a colum vector specifying different lengths
% * linkMode: can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)
% * input.rafLength: the length of the random access frame (RAF) (type: integer)
% * input.bitsPerSymbol: depends on the modulation scheme (type: integer)
% * input.fecRate: Forwar Error Correction rate; depends on the modulation scheme (type: double)
%
% Output
%
% output.queues: a logical matrix of input.sources rows and input.queueLength columns where each cell is set to 1 if the packet was successfully decoded or 0 if it was not
% output.delays: a matrix of input.sources rows and input.queueLength columns where each cell is set to the index of the RAF in which the packet was successfully decoded
% output.duration: the number of RAFs that have been generated to process all the input queues, needed to compute the average load and throughput

% TODO: complete the SUL mode first (1) [Issue: https://github.com/afcuttin/jsac/issues/29]
% TODO: complete the SDL mode second (2) [Issue: https://github.com/afcuttin/jsac/issues/31]
% TODO: complete the TDL mode third (3) [Issue: https://github.com/afcuttin/jsac/issues/30]
% TODO: complete the TUL mode fourth (4) [Issue: https://github.com/afcuttin/jsac/issues/28]

validateattributes(numberOfSources,{'numeric'},{'integer','positive'},mfilename,'numberOfSources',1);
validateattributes(queueLength,{'numeric'},{'vector','nonempty','integer','positive'},mfilename,'queueLength',2);
validatestring(linkMode,{'tul','sul','tdl','sdl'},mfilename,'linkMode',3);
assert(iscolumn(queueLength),'Variable queueLength must be a column vector');
assert((size(queueLength,1) == numberOfSources) || (size(queueLength,1) == 1),'The length of queueLength must be 1 or must be equal to numberOfSources');

input.sources             = numberOfSources;
input.queueLength         = queueLength;
input.linkMode            = linkMode;
% input.sinrThreshold     =
input.rafLength           = 7;
input.burstMaxRepetitions = 4;
input.bitsPerSymbol       = 3; % 8psk
input.fecRate             = 3/5;

queues.status        = input.queueLength .* ones(input.sources,1);
output.queues        = zeros(input.sources,max(input.queueLength));
output.delays        = zeros(input.sources,max(input.queueLength));
output.retries       = zeros(input.sources,max(input.queueLength));
outputMatrixSize     = size(output.queues);
% queste sono tutte le variabili ereditate dal vecchio codice, se possibile provvedere al refactoring
source.number        = input.sources;
% TODO: define a proper size of the RAF with respect to the number of actual sources [Issue: https://github.com/afcuttin/jsac/issues/4]
raf.length           = input.rafLength; % Casini et al., 2007, pag.1413
% simulationTime     = 1000; % total number of RAF
% simulationTime     = 10; % for TEST purposes only - comment when doing the real thing
% capturePar.status    = 3;
% sicPar.maxIter     = 16;
sicPar.maxIter       = 1;
sicPar.minIter       = 1;
% capturePar.criterion = 'power';
% capturePar.type      = 'basic';

% if strcmp(input.linkMode,'tul')
%     capturePar.threshold = input.bitsPerSymbol * input.fecRate;
% elseif strcmp(input.linkMode,'tdl') || strcmp(input.linkMode,'sdl') || strcmp(input.linkMode,'sul')
%     capturePar.threshold = 2^(input.bitsPerSymbol * input.fecRate) - 1;
% end

capturePar.threshold = 14; % the value of the parameter is obtained as follows: thr_val_dB + 11 = capturePar.threshold; thr_val_dB values are -10:1:10 % TEST: delete this line after testing

% the following for cycle should go away
for it = sicPar.minIter:sicPar.maxIter

    maxIter = it;

    % ackdPacketCount          = 0;
    % pcktTransmissionAttempts = 0;
    % pcktCollisionCount       = 0;
    source.status            = zeros(1,source.number);
    % legit source statuses are always non-negative integers and equal to:
    % 0: source has no packet ready to be transmitted (is idle)
    % 1: source has a packet ready to be transmitted, either because new data must be sent or a previously collided packet has waited the backoff time
    % integer greater than 1: source is backlogged due to previous packets collision, the integer value corresponds to the number of attempts made to get the latest burst acknowledged
    source.backoff           = zeros(1,source.number); % probably useless
    % pcktGenerationTimestamp  = zeros(1,source.number);
    output.duration          = 0;

        switch input.linkMode
            case 'tul' % random access method is Coded Slotted Aloha

                % carico il file che contiene le probabilità di cattura
                load('Captures_TUL_3','C_R_TUL_3','R_v','S_v');
                capturePar.rateThrVec      = R_v;
                capturePar.probability3seg = C_R_TUL_3;
                load('Captures_TUL_4','C_R_TUL_4');
                capturePar.probability4seg = C_R_TUL_4;
                capturePar.accessMethod    = 'csa';
                % TODO: update with correct capure probabilites after testing [Issue: https://github.com/afcuttin/jsac/issues/8]
                % load('capt_SUL.mat','C');
                % capturePar.probability = C;

                while sum(queues.status) > 0

                    assert(all(queues.status >= 0),'The status of a queue cannot be negative');

                    output.duration = output.duration + 1; % in multiples of RAF
                    raf.status      = zeros(source.number,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(source.number,raf.length);
                    % changedSlots    = 0;

                    % create the RAF
                    % TODO: update with correct number of bursts after testing [Issue: https://github.com/afcuttin/jsac/issues/2]
                    numberOfBursts = 3; % CSA method
                    % numberOfBursts = 2; % for testing purposes
                    for eachSource1 = 1:source.number
                        % TODO: inserire esperimento aleatorio per la scelta  del numero di pacchetti (come in IRSA) [Issue: https://github.com/afcuttin/jsac/issues/11]
                        % TODO: inserire la possibilità di fare arrivi di Poisson [Issue: https://github.com/afcuttin/jsac/issues/22]
                        if queues.status(eachSource1) > 0
                            if source.status(1,eachSource1) == 0 % a new burst can be sent
                                source.status(1,eachSource1)      = 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) >= 1 && source.status(1,eachSource1) <= input.burstMaxRepetitions  % backlogged source
                                source.status(1,eachSource1)      = source.status(1,eachSource1) + 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) > input.burstMaxRepetitions  % backlogged source, reached maximum number of attempts, discard backlogged burst
                                if queues.status(eachSource1) > 1
                                    queues.status(eachSource1)        = queues.status(eachSource1) - 1; % permanently drop unconfirmed packet
                                    % proceed with the transmission of the next packet in the queue
                                    source.status(1,eachSource1)      = 1;
                                    [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                    raf.status(eachSource1,pcktTwins) = 1;
                                    raf.twins(eachSource1,:)          = rafRow;
                                elseif queues.status(eachSource1) == 1
                                    % backlogged source, reached maximum number of attempts, discard backlogged burst, which is also the last in the queue
                                    queues.status(eachSource1)   = queues.status(eachSource1) - 1;
                                    assert(queues.status(eachSource1) >= 0,'negative queue status');
                                    source.status(1,eachSource1) = 0;
                                end
                            else
                                error('Unlegit sourse status.');
                            end
                        end
                    end

                    assert(all(sum(raf.status,2) == 3)) % TEST: delete this line after testing

                    acked.slot   = [];
                    acked.source = [];
                    iter         = 0;
                    enterTheLoop = true;

                    while (sum(raf.slotStatus) > 0 && iter <= maxIter) || enterTheLoop
                        enterTheLoop = false;
                        % decoding
                        [decRaf,decAcked] = decoding(raf,capturePar);
                        % update acked bursts list
                        acked.slot   = [acked.slot,decAcked.slot];
                        acked.source = [acked.source,decAcked.source];
                        % perform interference cancellation
                        icRaf = ic(decRaf,decAcked);
                        % start again
                        raf = icRaf;
                        iter = iter + 1;
                    end

                    % check for duplicates
                    count = histc(acked.source,unique(acked.source));
                    duplicatesExist = sum(count > 1) > 0;
                    assert(~duplicatesExist,'Error in the Successive Interference Cancellation process: one or more sources are acknowledged more than once');

                    % pcktTransmissionAttempts = pcktTransmissionAttempts + sum(source.status == 1); % "the normalized MAC load G does not take into account the replicas" Casini et al., 2007, pag.1411; "The performance parameter is throughput (measured in useful packets received per slot) vs. load (measured in useful packets transmitted per slot" Casini et al., 2007, pag.1415
                    % ackdPacketCount = ackdPacketCount + numel(acked.source);

                    % fprintf('Acked sources %f \n',acked.source); % TEST: delete this line after testing
                    % fprintf('Queues status %f \n',queues.status); % TEST: delete this line after testing
                    fprintf('Acked sources\n'); % TEST: delete this line after testing
                    acked.source % TEST: delete this line after testing
                    fprintf('Queues status\n'); % TEST: delete this line after testing
                    queues.status % TEST: delete this line after testing
                    % update the confirmed packets' status
                    % output.queues(sub2ind([input.sources max(input.queueLength)],transpose(acked.source),queues.status([acked.source]))) = 1;
                    output.queues(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = 1;
                    output.delays(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = output.duration;
                    % TODO: registrare il numero di ritrasmissioni necessarie per ogni pacchetto [Issue: https://github.com/afcuttin/jsac/issues/26]

                    % update the transmission queues
                    queues.status([acked.source]) = queues.status([acked.source]) - 1;
                    source.status([acked.source]) = 0; % update sources statuses
                    source.status(source.status < 0) = 0; % idle sources stay idle (see permitted statuses above)
                    % memoryless process (no retransmission attempts)
                    % queues.status = queues.status - 1;
                    % source.status = source.status - 1; % update sources statuses
                end

            case 'sul' % random access method is CRDSA

                % carico il file che contiene le probabilità di cattura
                % load('Captures_SUL','C_SUL','S_v'); % TEST: uncomment after testing
                % capturePar.sinrThrVec  = S_v; % TEST: uncomment after testing
                % capturePar.probability = C_SUL; % TEST: uncomment after testing
                % capturePar.accessMethod    = 'crdsa'; % TEST: uncomment after testing
                load('Captures_SUL','C_SUL');
                capturePar.probability = C_SUL;
                numberOfBursts = 2;
                capturePar.accessMethod = 'crdsa';

                while sum(queues.status) > 0

                    assert(all(queues.status >= 0),'The status of a queue cannot be negative');

                    output.duration = output.duration + 1; % in multiples of RAF
                    raf.status      = zeros(source.number,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(source.number,raf.length);

                    % TODO: completare Satellite Uplink [Issue: https://github.com/afcuttin/jsac/issues/3]

                    % create the RAF
                    for eachSource1 = 1:source.number
                        if queues.status(eachSource1) > 0
                            if source.status(1,eachSource1) == 0 % a new burst can be sent
                                source.status(1,eachSource1)      = 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                                % TODO: record here the time (raf id) when the burst has been transmitted for the first time (1 of 2) [Issue: https://github.com/afcuttin/jsac/issues/34]
                            elseif source.status(1,eachSource1) >= 1 && source.status(1,eachSource1) < input.burstMaxRepetitions  % backlogged source
                                source.status(1,eachSource1)      = source.status(1,eachSource1) + 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) >= input.burstMaxRepetitions  % backlogged source, reached maximum number of attempts, discard backlogged burst
                                if queues.status(eachSource1) > 1
                                    queues.status(eachSource1)        = queues.status(eachSource1) - 1; % permanently drop unconfirmed packet
                                    % proceed with the transmission of the next packet in the queue
                                    source.status(1,eachSource1)      = 1;
                                    [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                    raf.status(eachSource1,pcktTwins) = 1;
                                    raf.twins(eachSource1,:)          = rafRow;
                                    % TODO: record here the time (raf id) when the burst has been transmitted for the first time (2 of 2) [Issue: https://github.com/afcuttin/jsac/issues/36]
                                elseif queues.status(eachSource1) == 1
                                    % backlogged source, reached maximum number of attempts, discard backlogged burst, which is also the last in the queue
                                    queues.status(eachSource1)   = queues.status(eachSource1) - 1;
                                    assert(queues.status(eachSource1) >= 0,'negative queue status');
                                    source.status(1,eachSource1) = 0;
                                end
                            else
                                error('Unlegit sourse status.');
                            end
                        end
                    end

                    acked.slot   = [];
                    acked.source = [];
                    iter         = 0;
                    enterTheLoop = true;

                    while (sum(raf.slotStatus) > 0 && iter <= maxIter) || enterTheLoop
                        enterTheLoop = false;
                        % decoding
                        [decRaf,decAcked] = decoding(raf,capturePar);
                        % update acked bursts list
                        acked.slot   = [acked.slot,decAcked.slot];
                        acked.source = [acked.source,decAcked.source];
                        % perform interference cancellation
                        icRaf = ic(decRaf,decAcked);
                        % start again
                        raf = icRaf;
                        iter = iter + 1;
                    end

                    % check for duplicates
                    count = histc(acked.source,unique(acked.source));
                    duplicatesExist = sum(count > 1) > 0;
                    assert(~duplicatesExist,'Error in the Successive Interference Cancellation process: one or more sources are acknowledged more than once');

                    % pcktTransmissionAttempts = pcktTransmissionAttempts + sum(source.status == 1); % "the normalized MAC load G does not take into account the replicas" Casini et al., 2007, pag.1411; "The performance parameter is throughput (measured in useful packets received per slot) vs. load (measured in useful packets transmitted per slot" Casini et al., 2007, pag.1415
                    % ackdPacketCount = ackdPacketCount + numel(acked.source);

                    % fprintf('Acked sources %f \n',acked.source); % TEST: delete this line after testing
                    % fprintf('Queues status %f \n',queues.status); % TEST: delete this line after testing
                    fprintf('Acked sources\n'); % TEST: delete this line after testing
                    acked.source % TEST: delete this line after testing
                    fprintf('Queues status\n'); % TEST: delete this line after testing
                    queues.status % TEST: delete this line after testing
                    % update the confirmed packets' status
                    % output.queues(sub2ind([input.sources max(input.queueLength)],transpose(acked.source),queues.status([acked.source]))) = 1; % TEST: delete this line after testing
                    output.queues(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = 1;
                    output.delays(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = output.duration;
                    % TODO: to evaluate delay, the slot when the burst has been acked shall be recorded [Issue: https://github.com/afcuttin/jsac/issues/35]
                    for ii = 1:numel(acked.source)
                        output.delaySlot(acked.source(ii),queues.status(acked.source(ii))) = acked.slot(ii);
                    end
                    % TODO: registrare il numero di ritrasmissioni necessarie per ogni pacchetto [Issue: https://github.com/afcuttin/jsac/issues/26]
                    for ii = 1:numel(acked.source)
                        output.retries(acked.source(ii),queues.status(acked.source(ii))) = source.status(acked.source(ii));
                    end

                    % update the transmission queues
                    queues.status([acked.source]) = queues.status([acked.source]) - 1;
                    source.status % TEST: delete this line after testing
                    source.status([acked.source]) = 0; % update sources statuses
                    assert(all(source.status >= 0) && all(source.status <= input.burstMaxRepetitions)) % TEST: delete this line after testing
                    source.status(source.status < 0) = 0; % idle sources stay idle (see permitted statuses above)
                    % memoryless process (no retransmission attempts)
                    % queues.status = queues.status - 1;
                    % source.status = source.status - 1; % update sources statuses
                end

            case {'sdl','tdl'} % no random access, just capture threshold

                % TODO: completare Satellite Downlink [Issue: https://github.com/afcuttin/jsac/issues/6]
                % TODO: completare Terrestrial Downlink [Issue: https://github.com/afcuttin/jsac/issues/5]

                if strcmp(input.linkMode,'sdl')
                    load('Captures_SDL','C_SDL');
                    capturePar.probability = C_SDL;
                    capturePar.sinrThrVec  = S_v;
                elseif strcmp(input.linkMode,'tdl')
                    load('Captures_TDL','C_TDL');
                    capturePar.probability = C_TDL;
                    capturePar.sinrThrVec  = S_v;
                end

            otherwise
                error('Please select one of the availables link modes (tul, sul, sdl, tdl).');
        end

% FIXME: TODO: le code in uscita dei pacchetti nelle matrici sono disallineate (prova) [Issue: https://github.com/afcuttin/jsac/issues/27]
output.queues = fliplr(output.queues); % because packets are updated in reverse order
output.delays = fliplr(output.delays); % because packets are updated in reverse order

end

output.rafLength = input.rafLength;
