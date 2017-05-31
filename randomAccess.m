function [outQueues,outDelays,outRetries,outFirstTx,outDuration,outRafLength,output] = randomAccess(numberOfSources,queueLength,linkMode,inputPars)
% function [outQueues,outDelays,outRetries,outFirstTx,outDuration,outRafLength,output] = randomAccess(numberOfSources,queueLength,linkMode)
%
% Simulation of Multiple Random Access
%
% *** Input parameters:
%
% * numberOfSources: number of sources (type: integer)
% * queueLength:     number of packets (burst) in the fifo queue for each source. The length must be the same for all sources (type: integer). Can be a scalar (all the queues have the same length) or a colum vector specifying different lengths
% * linkMode:        can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
%
% The inputPars structure is not mandatory. If not specified, default parameters will be used.
%
% * inputPars.sinrThreshold:    value of the SINR threshold to be used (type: integer, default 4)
% * inputPars.rafLength:        the length of the random access frame (RAF) (type: integer)
% * inputPars.poissonThreshold: threshold for a random experiment that enables sources to transimt (type: real, default 1)
% * inputPars.retryLimit:       number of transmission retries for a given data packet (type: integer, default 4)
% * inputPars.timeConstraint:   the maximum allowed duration of a simulation, measured in slots (type: integer, default inf)
%
% *** Outputs:
%
% * outQueues:    a logical matrix of numberOfSources rows and queueLength columns where each cell is set to 1 if the packet was successfully decoded or 0 if it was not
% * outDelays:    a matrix of numberOfSources rows and queueLength columns where each cell is set to the index of the RAF in which the packet was successfully decoded
% * outRetries:   a matrix of numberOfSources rows and queueLength columns where each cell is set to the number of times a packet has been sent in a RAF
% * outFirstTx:   a matrix of numberOfSources rows and queueLength columns where each cell is set to the index of the RAF when the data packet has been transmitted for the first time
% * outDuration:  the number of RAFs that have been generated to process all the input queues, needed to compute the average load and throughput
% * outRafLength: the number of slots present in a RAF, used to compute the duration in slot

% TODO: assess the need of code deduplication (6) [Issue: https://github.com/afcuttin/jsac/issues/52]

validateattributes(numberOfSources,{'numeric'},{'integer','positive'},mfilename,'numberOfSources',1);
validateattributes(queueLength,{'numeric'},{'vector','nonempty','integer','positive'},mfilename,'queueLength',2);
validatestring(linkMode,{'tul','sul','tdl','sdl'},mfilename,'linkMode',3);
assert(iscolumn(queueLength),'Variable queueLength must be a column vector');
assert((size(queueLength,1) == numberOfSources) || (size(queueLength,1) == 1),'The length of queueLength must be 1 or must be equal to numberOfSources');

if exist('inputPars','var')
    if ~isfield(inputPars,'poissonThreshold')
        inputPars.poissonThreshold = 1;
    elseif isfield(inputPars,'poissonThreshold')
        validateattributes(inputPars.poissonThreshold,{'numeric'},{'nonempty','nonzero','>',0,'<=',1},mfilename,'inputPars.poissonThreshold',4);
    end
    if isfield(inputPars,'rafLength')
        raf.length = inputPars.rafLength;
    elseif ~isfield(inputPars,'rafLength')
        raf.length = 10;
    end
    if ~isfield(inputPars,'retryLimit')
        inputPars.retryLimit = 4;
    end
    if ~isfield(inputPars,'sinrThreshold')
        inputPars.sinrThreshold = 4;
    end
    if ~isfield(inputPars,'timeConstraint')
        inputPars.timeConstraint = inf(1);
    end
elseif ~exist('inputPars','var')
    inputPars.poissonThreshold = 1;
    inputPars.retryLimit       = 4;
    inputPars.sinrThreshold    = 4;
% TODO: define a proper size of the RAF with respect to the number of actual sources [Issue: https://github.com/afcuttin/jsac/issues/4]
    raf.length                 = 10;
    inputPars.timeConstraint   = inf(1);
end
% queste sono tutte le variabili ereditate dal vecchio codice, se possibile provvedere al refactoring
sicPar.maxIter       = 2;
sicPar.minIter       = 1;

queueLength      = queueLength .* ones(numberOfSources,1); % in any case, queueLength becomes a column vector
queues.status    = ones(numberOfSources,1);
output.queues    = zeros(numberOfSources,max(queueLength));
output.delays    = zeros(numberOfSources,max(queueLength));
output.delaySlot = zeros(numberOfSources,max(queueLength));
output.retries   = zeros(numberOfSources,max(queueLength));
output.firstTx   = zeros(numberOfSources,max(queueLength));
output.duration  = 0;
outputMatrixSize = size(output.queues);
source.status    = zeros(1,numberOfSources);
% legit source statuses are always non-negative integers and equal to:
% 0: source has no packet ready to be transmitted (is idle)
% 1: source has a packet ready to be transmitted, either because new data must be sent or a previously collided packet has waited the backoff time
% integer greater than 1: source is backlogged due to previous packets collision, the integer value corresponds to the number of attempts made to get the latest burst acknowledged

        switch linkMode
            case 'tul' % random access method is Coded Slotted Aloha
                % carico il file che contiene le probabilità di cattura
                load('Captures_TUL_3','C_R_TUL_3','R_v');
                capturePar.rateThrVec      = R_v;
                capturePar.probability3seg = C_R_TUL_3;
                load('Captures_TUL_4','C_R_TUL_4');
                capturePar.probability4seg = C_R_TUL_4;
                capturePar.accessMethod    = 'csa'; % NOTE: this setting can't be controlled from the output

                while any(queues.status <= queueLength) && output.duration * raf.length < inputPars.timeConstraint

                    assert(all(queues.status <= queueLength+1),'The number of confirmed packets shall not exceed the lenght of the queue.');

                    output.duration = output.duration + 1; % in multiples of RAF
                    raf.status      = zeros(numberOfSources,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(numberOfSources,raf.length);

                    % create the RAF
                    % NOTE: prima di partire col ciclo, trovare le sorgenti che hanno ancora pacchetti in coda da smaltire, e ciclare solo su quelle, così si può eliminare il condizionale di 116 (4 righe più in basso)
                    % NOTE: the following for cycle is used here and in the other link mode: it can be converted in a single function to prevent code duplication and problems in its update
                    if ~isfield(inputPars,'poissonThreshold') || inputPars.poissonThreshold == 1
                        enabledSources = ones(numberOfSources,1);
                    elseif isfield(inputPars,'poissonThreshold')
                        enabledSources = rand(numberOfSources,1) <= inputPars.poissonThreshold;
                    end
                    for eachSource1 = 1:numberOfSources
                        pcktRepExp = rand(1);
                        if pcktRepExp <= 1/3
                            numberOfBursts = 4;
                        elseif pcktRepExp <= (1/3 + 2/3)
                            numberOfBursts = 3;
                        end
                        if queues.status(eachSource1) <= queueLength(eachSource1)
                            if source.status(1,eachSource1) == 0 && enabledSources(eachSource1) == 1 % a new burst can be sent
                                source.status(1,eachSource1)      = 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                                output.firstTx(eachSource1,queues.status(eachSource1)) = output.duration;
                            elseif source.status(1,eachSource1) == 0 && enabledSources(eachSource1) == 0
                                % stay idle
                                % NOTE: the conditional on the probability is introduced in order to have a Poisson distribution of the arrivals (that is, transmission). However, the approach followed here does not result in a pure Poisson distribution. In fact, the random experiment is ran only for those sources that are idle (that is: they have a new packet to be sent) or have exceeded the retry limit. Backlogged sources are not subject to the random experiment and therefore can transmit their packet until it is acknowledged or the retry limit is exceeded. This is done to have a behaviour that is as close as possible to the way a real device works.
                            elseif source.status(1,eachSource1) >= 1 && source.status(1,eachSource1) < inputPars.retryLimit  % backlogged source
                                source.status(1,eachSource1)      = source.status(1,eachSource1) + 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) >= inputPars.retryLimit  % backlogged source, reached maximum retry limit, discard backlogged burst
                                if queues.status(eachSource1) < queueLength(eachSource1)
                                    queues.status(eachSource1)        = queues.status(eachSource1) + 1; % permanently drop unconfirmed packet
                                    % proceed with the transmission of the next packet in the queue
                                    if enabledSources(eachSource1) == 1
                                        source.status(1,eachSource1)      = 1;
                                        [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                        raf.status(eachSource1,pcktTwins) = 1;
                                        raf.twins(eachSource1,:)          = rafRow;
                                        output.firstTx(eachSource1,queues.status(eachSource1)) = output.duration;
                                    elseif enabledSources(eachSource1) == 0
                                        source.status(1,eachSource1) = 0;
                                    end
                                elseif queues.status(eachSource1) == queueLength(eachSource1)
                                    % backlogged source, reached maximum number of attempts, discard the backlogged burst, which is also the last in the queue
                                    queues.status(eachSource1)   = queues.status(eachSource1) + 1;
                                    source.status(1,eachSource1) = 0;
                                end
                            else
                                error('Unlegit source status.');
                            end
                        end
                    end

                    acked.slot   = [];
                    acked.source = [];

                    switch capturePar.accessMethod
                        case 'csa-p'
                            % decoding
                            [decRaf,decAcked] = decoding(raf,capturePar);
                            % update acked bursts list
                            acked.slot   = [acked.slot,decAcked.slot];
                            acked.source = [acked.source,decAcked.source];
                        case 'csa-pip'
                            % decoding
                            [decRaf,decAcked] = decoding(raf,capturePar);
                            % update acked bursts list
                            acked.slot   = [acked.slot,decAcked.slot];
                            acked.source = [acked.source,decAcked.source];
                            % perform interference cancellation (only once)
                            icRaf = ic(decRaf,decAcked);
                            % start again
                            raf = icRaf;
                            % decoding
                            [decRaf,decAcked] = decoding(raf,capturePar);
                            % update acked bursts list
                            acked.slot   = [acked.slot,decAcked.slot];
                            acked.source = [acked.source,decAcked.source];
                        case 'csa'
                            iter         = 0;
                            enterTheLoop = true;
                            while (sum(raf.slotStatus) > 0 && iter <= sicPar.maxIter) || enterTheLoop
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
                        otherwise
                            error('Please select one of the availables CSA modes.');
                    end

                    % check for duplicates
                    count = histc(acked.source,unique(acked.source));
                    duplicatesExist = sum(count > 1) > 0;
                    assert(~duplicatesExist,'Error in the Successive Interference Cancellation process: one or more sources are acknowledged more than once');

                    % update the confirmed packets' status
                    output.queues(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = 1;
                    output.delays(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = output.duration;
                    for ii = 1:numel(acked.source)
                        output.retries(acked.source(ii),queues.status(acked.source(ii))) = source.status(acked.source(ii));
                    end
                    for ii = 1:numel(acked.source)
                        output.delaySlot(acked.source(ii),queues.status(acked.source(ii))) = acked.slot(ii);
                    end

                    % update the transmission queues
                    queues.status([acked.source]) = queues.status([acked.source]) + 1;
                    source.status([acked.source]) = 0; % update sources statuses
                    source.status(source.status < 0) = 0; % idle sources stay idle (see permitted statuses above)
                    % memoryless process (no retransmission attempts) NOTE: this is probably equivalent to setting inputPars.retryLimit = 1
                    % queues.status = queues.status - 1;
                    % source.status = source.status - 1; % update sources statuses
                end
                output.rafLength = raf.length;

            case 'sul' % random access method is CRDSA

                % carico il file che contiene le probabilità di cattura
                load('Captures_SUL','C_SUL','S_v');
                capturePar.probability    = C_SUL;
                [~,capturePar.sinrThrInd] = min(abs(S_v - inputPars.sinrThreshold));
                capturePar.accessMethod   = 'crdsa';
                numberOfBursts            = 2;

                while any(queues.status <= queueLength) && output.duration * raf.length < inputPars.timeConstraint

                    assert(all(queues.status <= queueLength+1),'The number of confirmed packets shall not exceed the lenght of the queue.');

                    output.duration = output.duration + 1; % in multiples of RAF
                    raf.status      = zeros(numberOfSources,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(numberOfSources,raf.length);

                    % create the RAF
                    % NOTE: prima di partire col ciclo, trovare le sorgenti che hanno ancora pacchetti in coda da smaltire, e ciclare solo su quelle, così si può eliminare il condizionale di 116 (4 righe più in basso)
                    if ~isfield(inputPars,'poissonThreshold') || inputPars.poissonThreshold == 1
                        enabledSources = ones(numberOfSources,1);
                    elseif isfield(inputPars,'poissonThreshold')
                        enabledSources = rand(numberOfSources,1) <= inputPars.poissonThreshold;
                    end
                    for eachSource1 = 1:numberOfSources
                        if queues.status(eachSource1) <= queueLength(eachSource1)
                            if source.status(1,eachSource1) == 0 && enabledSources(eachSource1) == 1 % a new burst can be sent
                                source.status(1,eachSource1)      = 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                                output.firstTx(eachSource1,queues.status(eachSource1)) = output.duration;
                            elseif source.status(1,eachSource1) == 0 && enabledSources(eachSource1) == 0
                                % stay idle
                                % NOTE: the conditional on the probability is introduced in order to have a Poisson distribution of the arrivals (that is, transmission). However, the approach followed here does not result in a pure Poisson distribution. In fact, the random experiment is ran only for those sources that are idle (that is: they have a new packet to be sent) or have exceeded the retry limit. Backlogged sources are not subject to the random experiment and therefore can transmit their packet until it is acknowledged or the retry limit is exceeded. This is done to have a behaviour that is as close as possible to the way a real device works.
                            elseif source.status(1,eachSource1) >= 1 && source.status(1,eachSource1) < inputPars.retryLimit  % backlogged source
                                source.status(1,eachSource1)      = source.status(1,eachSource1) + 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) >= inputPars.retryLimit  % backlogged source, reached maximum retry limit, discard backlogged burst
                                if queues.status(eachSource1) < queueLength(eachSource1)
                                    queues.status(eachSource1)        = queues.status(eachSource1) + 1; % permanently drop unconfirmed packet
                                    % proceed with the transmission of the next packet in the queue
                                    if enabledSources(eachSource1) == 1
                                        source.status(1,eachSource1)      = 1;
                                        [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                        raf.status(eachSource1,pcktTwins) = 1;
                                        raf.twins(eachSource1,:)          = rafRow;
                                        output.firstTx(eachSource1,queues.status(eachSource1)) = output.duration;
                                    elseif enabledSources(eachSource1) == 0
                                        source.status(1,eachSource1) = 0;
                                    end
                                elseif queues.status(eachSource1) == queueLength(eachSource1)
                                    % backlogged source, reached maximum number of attempts, discard backlogged burst, which is also the last in the queue
                                    queues.status(eachSource1)   = queues.status(eachSource1) + 1;
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

                    while (sum(raf.slotStatus) > 0 && iter <= sicPar.maxIter) || enterTheLoop
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

                    % update the confirmed packets' status
                    output.queues(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = 1;
                    output.delays(sub2ind(outputMatrixSize,transpose(acked.source),queues.status([acked.source]))) = output.duration;
                    for ii = 1:numel(acked.source)
                        output.retries(acked.source(ii),queues.status(acked.source(ii))) = source.status(acked.source(ii));
                    end
                    for ii = 1:numel(acked.source)
                        output.delaySlot(acked.source(ii),queues.status(acked.source(ii))) = acked.slot(ii);
                    end

                    % update the transmission queues
                    queues.status([acked.source]) = queues.status([acked.source]) + 1;
                    source.status([acked.source]) = 0; % update sources statuses
                    assert(all(source.status >= 0) && all(source.status <= inputPars.retryLimit)) % NOTE: this check on the statuses is a duplicate of the one performed above: "error('Unlegit sourse status.')"
                    source.status(source.status < 0) = 0; % idle sources stay idle (see permitted statuses above)
                    % memoryless process (no retransmission attempts) NOTE: this is probably equivalent to setting inputPars.retryLimit = 1
                    % queues.status = queues.status + 1;
                    % source.status = source.status - 1; % update sources statuses
                end
                output.rafLength = raf.length;

            case {'sdl','tdl'} % no random access, just capture threshold

                switch linkMode
                    case 'sdl'
                        load('Captures_SDL');
                        capturePar.probability = C_SDL;
                        capturePar.sinrThrVec  = S_v;
                    case 'tdl'
                        load('Captures_TDL');
                        capturePar.probability = C_TDL;
                        capturePar.sinrThrVec  = S_v;
                    otherwise
                        % no need to throw an error
                end
                [~,sinrThrInd] = min(abs(capturePar.sinrThrVec - inputPars.sinrThreshold));

                while any(queues.status <= queueLength) && output.duration * raf.length < inputPars.timeConstraint

                    assert(all(queues.status <= queueLength+1),'The number of confirmed packets shall not exceed the lenght of the queue.');
                    output.duration = output.duration + 1; % in multiples of RAF

                    % find idle, backlogged and unsuccessful sources
                    inactiveSources     = find(queues.status > queueLength);
                    idleSources         = find(source.status == 0);
                    idleSources         = setdiff(idleSources,inactiveSources);
                    backloggedSources   = find(ismember(source.status,[1:1:inputPars.retryLimit]));
                    backloggedSources   = setdiff(backloggedSources,inactiveSources);
                    unsuccessfulSources = find(source.status == inputPars.retryLimit + 1);
                    atEndOfQueue        = find(queues.status == queueLength);
                    atEndOfQueue        = intersect(unsuccessfulSources,atEndOfQueue);
                    unsuccessfulSources = setdiff(unsuccessfulSources,[inactiveSources ; atEndOfQueue]);
                    assert(all(source.status <= inputPars.retryLimit + 1),'A source status is one unit too big');

                    % TODO: if inputPars.poissonThreshold is enabled, sdl and tdl fail the test (1) [Issue: https://github.com/afcuttin/jsac/issues/53]
                    if ~isfield(inputPars,'poissonThreshold') || inputPars.poissonThreshold == 1
                        enabledSources = idleSources;
                    elseif isfield(inputPars,'poissonThreshold')
                        % enabledSources = rand(numberOfSources,1) <= inputPars.poissonThreshold; % così è nei casi precedenti
                        enabledSources = idleSources(find([rand(numel(idleSources),1) <= inputPars.poissonThreshold] ));
                    end
                    % update the status of idle and unsuccessful sources
                    source.status(enabledSources)         = 1;
                    source.status(unsuccessfulSources) = 1; % unsuccessful sources drop the current packet and move to the next one
                    source.status(atEndOfQueue)        = 0; % unsuccessful sources at the end of the queue drop the current packet and stay permanently idle
                    queues.status(unsuccessfulSources) = queues.status(unsuccessfulSources) + 1;
                    queues.status(atEndOfQueue)        = queues.status(atEndOfQueue) + 1;

                    % update the firstTx matrix
                    if ~isempty(enabledSources)
                        output.firstTx(sub2ind(outputMatrixSize,transpose(enabledSources),queues.status(enabledSources)))                 = output.duration;
                    end
                    if ~isempty(unsuccessfulSources)
                        output.firstTx(sub2ind(outputMatrixSize,transpose(unsuccessfulSources),queues.status(unsuccessfulSources))) = output.duration;
                    end

                    % run the random experiment to determine successful packet reception (acknowledgment)
                    randomExperiments      = rand(numberOfSources,1);
                    successfulTransmission = randomExperiments <= capturePar.probability(sinrThrInd);
                    ackedSources           = find(successfulTransmission == 1);
                    ackedSources           = setdiff(ackedSources,[inactiveSources ; atEndOfQueue]); % NOTE: it is probably better to setdiff between ackedSources and enabledSources
                    backlSources           = find(~successfulTransmission == 1);
                    backlSources           = setdiff(backlSources,[inactiveSources ; atEndOfQueue]); % NOTE: it is probably better to setdiff between ackedSources and enabledSources

                    % update output matrices and transmission queues for acked sources
                    if ~isempty(ackedSources)
                        output.queues(sub2ind(outputMatrixSize,ackedSources,queues.status(ackedSources))) = 1;
                        output.delays(sub2ind(outputMatrixSize,ackedSources,queues.status(ackedSources))) = output.duration;
                        for ii = 1:numel(ackedSources)
                            output.retries(ackedSources(ii),queues.status(ackedSources(ii))) = source.status(ackedSources(ii));
                        end
                    end

                    queues.status(ackedSources) = queues.status(ackedSources) + 1;
                    source.status(ackedSources) = 0; % update sources statuses
                    source.status(backlSources) = source.status(backlSources) + 1; % update sources statuses
                end
                output.rafLength = 1;
            otherwise
                error('Please select one of the availables link modes (tul, sul, sdl, tdl).');
        end
outQueues    = output.queues;
outDelays    = output.delays;
outRetries   = output.retries;
outFirstTx   = output.firstTx;
outDuration  = output.duration;
outRafLength = output.rafLength;
