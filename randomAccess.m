function [outQueues,outDelays,outRetries,outFirstTx,outDuration,outRafLength,output] = randomAccess(numberOfSources,queueLength,linkMode)
% function [outQueues,outDelays,outRetries,outFirstTx,outDuration,outRafLength,output] = randomAccess(numberOfSources,queueLength,linkMode)
%
% Simulation of Multiple Random Access
%
% Input parameters:
% * numberOfSources: number of sources (type: integer)
% * queueLength: number of packets (burst) in the fifo queue for each source. The length must be the same for all sources (type: integer). Can be a scalar (all the queues have the same length) or a colum vector specifying different lengths
% * linkMode: can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)
% * input.rafLength: the length of the random access frame (RAF) (type: integer) NOTE: this input parameter is deprecated. It will be specified in another way
% * input.bitsPerSymbol: depends on the modulation scheme (type: integer)
% * input.fecRate: Forwar Error Correction rate; depends on the modulation scheme (type: double)
% TODO: update randomAccess.m function help with correct type of inputs (99) [Issue: https://github.com/afcuttin/jsac/issues/51]
%
% Output
% * outQueues: a logical matrix of numberOfSources rows and queueLength columns where each cell is set to 1 if the packet was successfully decoded or 0 if it was not
% * outDelays: a matrix of numberOfSources rows and queueLength columns where each cell is set to the index of the RAF in which the packet was successfully decoded
% * outRetries:
% * outFirstTx:
% * outDuration:
% * outRafLength: the number of RAFs that have been generated to process all the input queues, needed to compute the average load and throughput

% TODO: assess the need of code deduplication (6) [Issue: https://github.com/afcuttin/jsac/issues/52]

validateattributes(numberOfSources,{'numeric'},{'integer','positive'},mfilename,'numberOfSources',1);
validateattributes(queueLength,{'numeric'},{'vector','nonempty','integer','positive'},mfilename,'queueLength',2);
validatestring(linkMode,{'tul','sul','tdl','sdl'},mfilename,'linkMode',3);
assert(iscolumn(queueLength),'Variable queueLength must be a column vector');
assert((size(queueLength,1) == numberOfSources) || (size(queueLength,1) == 1),'The length of queueLength must be 1 or must be equal to numberOfSources');

% NOTE: it could be good that the following input parameter could be configured from outside the funcion (maybe with a different configuration script that generates some .m file from which the variables are loaded)
input.sources             = numberOfSources;
input.linkMode            = linkMode;
input.sinrThreshold       = 4; % value in dB NOTE: this parameter is no longer used
input.burstMaxRepetitions = 4; % NOTE: this is the retry limit, maybe rename to input.retryLimit
input.bitsPerSymbol       = 3; % 8psk % NOTE: this parameter is no longer used
input.fecRate             = 3/5; % NOTE: this parameter is no longer used

queueLength          = queueLength .* ones(input.sources,1); % in any case, queueLength becomes a column vector
queues.status        = ones(input.sources,1);
output.queues        = zeros(input.sources,max(queueLength));
output.delays        = zeros(input.sources,max(queueLength));
output.retries       = zeros(input.sources,max(queueLength));
output.firstTx       = zeros(input.sources,max(queueLength));
outputMatrixSize     = size(output.queues);
% queste sono tutte le variabili ereditate dal vecchio codice, se possibile provvedere al refactoring
source.number        = input.sources;
% TODO: define a proper size of the RAF with respect to the number of actual sources [Issue: https://github.com/afcuttin/jsac/issues/4]
raf.length           = 10;
sicPar.maxIter       = 1;
sicPar.minIter       = 1;

source.status = zeros(1,source.number);
% legit source statuses are always non-negative integers and equal to:
% 0: source has no packet ready to be transmitted (is idle)
% 1: source has a packet ready to be transmitted, either because new data must be sent or a previously collided packet has waited the backoff time
% integer greater than 1: source is backlogged due to previous packets collision, the integer value corresponds to the number of attempts made to get the latest burst acknowledged
output.duration = 0;

        switch input.linkMode
            case 'tul' % random access method is Coded Slotted Aloha
                % carico il file che contiene le probabilità di cattura
                load('Captures_TUL_3','C_R_TUL_3','R_v');
                capturePar.rateThrVec      = R_v;
                capturePar.probability3seg = C_R_TUL_3;
                load('Captures_TUL_4','C_R_TUL_4');
                capturePar.probability4seg = C_R_TUL_4;
                capturePar.accessMethod    = 'csa'; % NOTE: this setting can't be controlled from the output

                while any(queues.status <= queueLength)

                    assert(all(queues.status <= queueLength+1),'The number of confirmed packets shall not exceed the lenght of the queue.');

                    output.duration = output.duration + 1; % in multiples of RAF
                    raf.status      = zeros(source.number,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(source.number,raf.length);

                    % create the RAF
                    % NOTE: prima di partire col ciclo, trovare le sorgenti che hanno ancora pacchetti in coda da smaltire, e ciclare solo su quelle, così si può eliminare il condizionale di 116 (4 righe più in basso)
                    % NOTE: the following for cycle is used here and in the other link mode: it can be converted in a single function to prevent code duplication and problems in its update
                    if ~exist('sogliaPoisson','var') || sogliaPoisson == 1
                        enabledSources = ones(source.number,1);
                    elseif exist('sogliaPoisson','var')
                        enabledSources = rand(source.number,1) <= sogliaPoisson;
                    end
                    for eachSource1 = 1:source.number
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
                            elseif source.status(1,eachSource1) >= 1 && source.status(1,eachSource1) < input.burstMaxRepetitions  % backlogged source
                                source.status(1,eachSource1)      = source.status(1,eachSource1) + 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) >= input.burstMaxRepetitions  % backlogged source, reached maximum retry limit, discard backlogged burst
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
                    % memoryless process (no retransmission attempts) NOTE: this is probably equivalent to setting input.burstMaxRepetitions = 1
                    % queues.status = queues.status - 1;
                    % source.status = source.status - 1; % update sources statuses
                end
                output.rafLength = raf.length;

            case 'sul' % random access method is CRDSA

                % carico il file che contiene le probabilità di cattura
                load('Captures_SUL','C_SUL','S_v');
                capturePar.probability    = C_SUL;
                [~,capturePar.sinrThrInd] = min(abs(S_v - input.sinrThreshold));
                capturePar.accessMethod   = 'crdsa';
                numberOfBursts            = 2;

                while any(queues.status <= queueLength)

                    assert(all(queues.status <= queueLength+1),'The number of confirmed packets shall not exceed the lenght of the queue.');

                    output.duration = output.duration + 1; % in multiples of RAF
                    raf.status      = zeros(source.number,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(source.number,raf.length);

                    % create the RAF
                    % NOTE: prima di partire col ciclo, trovare le sorgenti che hanno ancora pacchetti in coda da smaltire, e ciclare solo su quelle, così si può eliminare il condizionale di 116 (4 righe più in basso)
                    if ~exist('sogliaPoisson','var') || sogliaPoisson == 1
                        enabledSources = ones(source.number,1);
                    elseif exist('sogliaPoisson','var')
                        enabledSources = rand(source.number,1) <= sogliaPoisson;
                    end
                    for eachSource1 = 1:source.number
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
                            elseif source.status(1,eachSource1) >= 1 && source.status(1,eachSource1) < input.burstMaxRepetitions  % backlogged source
                                source.status(1,eachSource1)      = source.status(1,eachSource1) + 1;
                                [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                                raf.status(eachSource1,pcktTwins) = 1;
                                raf.twins(eachSource1,:)          = rafRow;
                            elseif source.status(1,eachSource1) >= input.burstMaxRepetitions  % backlogged source, reached maximum retry limit, discard backlogged burst
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
                    assert(all(source.status >= 0) && all(source.status <= input.burstMaxRepetitions)) % NOTE: this check on the statuses is a duplicate of the one performed above: "error('Unlegit sourse status.')"
                    source.status(source.status < 0) = 0; % idle sources stay idle (see permitted statuses above)
                    % memoryless process (no retransmission attempts) NOTE: this is probably equivalent to setting input.burstMaxRepetitions = 1
                    % queues.status = queues.status + 1;
                    % source.status = source.status - 1; % update sources statuses
                end
                output.rafLength = raf.length;

            case {'sdl','tdl'} % no random access, just capture threshold

                switch input.linkMode
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
                [~,sinrThrInd] = min(abs(capturePar.sinrThrVec - input.sinrThreshold));

                while any(queues.status <= queueLength)

                    assert(all(queues.status <= queueLength+1),'The number of confirmed packets shall not exceed the lenght of the queue.');
                    output.duration = output.duration + 1; % in multiples of RAF

                    % find idle, backlogged and unsuccessful sources
                    inactiveSources     = find(queues.status > queueLength);
                    idleSources         = find(source.status == 0);
                    idleSources         = setdiff(idleSources,inactiveSources);
                    backloggedSources   = find(ismember(source.status,[1:1:input.burstMaxRepetitions]));
                    backloggedSources   = setdiff(backloggedSources,inactiveSources);
                    unsuccessfulSources = find(source.status == input.burstMaxRepetitions + 1);
                    atEndOfQueue        = find(queues.status == queueLength);
                    atEndOfQueue        = intersect(unsuccessfulSources,atEndOfQueue);
                    unsuccessfulSources = setdiff(unsuccessfulSources,[inactiveSources ; atEndOfQueue]);
                    assert(all(source.status <= input.burstMaxRepetitions + 1),'A source status is one unit too big');

                    if ~exist('sogliaPoisson','var') || sogliaPoisson == 1
                        enabledSources = idleSources;
                    elseif exist('sogliaPoisson','var')
                        enabledSources = idleSources(find([rand(numel(idleSources),1) <= sogliaPoisson] ))
                    end
                    % update the status of idle and unsuccessful sources
                    source.status(enabledSources)         = 1;
                    source.status(unsuccessfulSources) = 1; % unsuccessful sources drop the current packet and move to the next one
                    source.status(atEndOfQueue)        = 0; % unsuccessful sources at the end of the queue drop the current packet and stay permanently idle
                    queues.status(unsuccessfulSources) = queues.status(unsuccessfulSources) + 1;
                    queues.status(atEndOfQueue)        = queues.status(atEndOfQueue) + 1;

                    % update the firstTx matrix
                    output.firstTx(sub2ind(outputMatrixSize,transpose(enabledSources),queues.status(enabledSources)))                 = output.duration;
                    output.firstTx(sub2ind(outputMatrixSize,transpose(unsuccessfulSources),queues.status(unsuccessfulSources))) = output.duration;

                    % run the random experiment to determine successful packet reception (acknowledgment)
                    randomExperiments      = rand(input.sources,1);
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
