function [output] = randomAccess(input)
% function [in] = randomAccess(out)
%
% Simulation of Multiple Random Access
%
% Input parameters:
% * input.sources: number of sources (type: integer)
% * input.queueLength: number of packets (burst) in the fifo queue for each source. The length must be the same for all sources (type: integer)
% * input.linkMode: can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)
% * input.rafLength: the length of the random access frame (RAF) (type: integer)
%
% Output
%
% output.queues: a logical matrix of input.sources rows and input.queueLength columns where each cell is set to 1 if the packet was successfully decoded or 0 if it was not
% output.duration: the number of RAFs that have been generated to process all the input queues

queues.status        = input.queueLength * ones(input.sources,1);
output.matrix        = zeros(input.sources,input.queueLength);
% queste sono tutte le variabili ereditate dal vecchio codice, se possibile provvedere al refactoring
source.number        = input.sources;
% TODO: define a proper size of the RAF with respect to the number of actual sources [Issue: https://github.com/afcuttin/jsac/issues/4]
raf.length           = input.rafLength; % Casini et al., 2007, pag.1413
% simulationTime     = 1000; % total number of RAF
% simulationTime     = 10; % for TEST purposes only - comment when doing the real thing
% capturePar.status    = 3;
% sicPar.maxIter     = 16;
sicPar.maxIter       = 4;
sicPar.minIter       = 4;
% capturePar.criterion = 'power';
% capturePar.type      = 'basic';

% TODO: define modulation key parameters (bitrate, bandwidth, modulation scheme, code rate) [Issue: https://github.com/afcuttin/jsac/issues/9]
% TODO: introduce the evaluation of the capture threshold as a function of the code rate and the number of bits per symbol of the selected modulation [Issue: https://github.com/afcuttin/jsac/issues/7]
capturePar.threshold = 14; % the value of the parameter is obtained as follows: thr_val_dB + 11 = capturePar.threshold; thr_val_dB values are -10:1:10

% the following for cycle should go away
for it = sicPar.minIter:sicPar.maxIter

    maxIter = it;

    ackdPacketCount          = 0;
    pcktTransmissionAttempts = 0;
    pcktCollisionCount       = 0;
    source.status            = zeros(1,source.number);
    % legit source statuses are always non-negative integers and equal to:
    % 0: source has no packet ready to be transmitted (is idle)
    % 1: source has a packet ready to be transmitted, either because new data must be sent or a previously collided packet has waited the backoff time
    % integer greater than 1: source is backlogged due to previous packets collision, the integer value corresponds to the number of attempts made to get the latest burst acknowledged
    source.backoff           = zeros(1,source.number); % probably useless
    pcktGenerationTimestamp  = zeros(1,source.number);
    output.duration          = 0;

        switch input.linkMode
            case 'tul' % random access method is Coded Slotted Aloha

                % carico il file che contiene le probabilità di cattura
                % load('Captures_TUL','C_TUL');
                % capturePar.probability = C_TUL;
                % TODO: update with correct capure probabilites after testing [Issue: https://github.com/afcuttin/jsac/issues/8]
                load('capt_SUL.mat','C');
                capturePar.probability = C;

                while sum(queues.status) > 0

                    assert(all(queues.status >= 0),'The status of a queue cannot be negative');

                    output.duration = output.duration + 1;
                    raf.status      = zeros(source.number,raf.length); % memoryless
                    raf.slotStatus  = int8(zeros(1,raf.length));
                    raf.twins       = cell(source.number,raf.length);
                    changedSlots    = 0;

                % create the RAF
                % TODO: update with correct number of bursts after testing [Issue: https://github.com/afcuttin/jsac/issues/2]
                numberOfBursts = 3; % CSA method
                numberOfBursts = 2; % for testing purposes
                for eachSource1 = 1:source.number
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
                            % decidere se fare qui lo scarto del pacchetto, o più avanti, quando si fa la verifica di quelli confermati
                            % marcare a 0 la posizione corrispondente al pacchetto scartato
                            % decrementare il contatore della coda dei pacchetti in ingresso
                            % procedere con una nuova trasmissione
                            source.status(1,eachSource1)      = 1;
                            [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                            raf.status(eachSource1,pcktTwins) = 1;
                            raf.twins(eachSource1,:)          = rafRow;
                        end
                    end
                end

                % acked.slot   = [];
                acked.source = [];
                iter         = 0;
                enterTheLoop = true;

                while (sum(raf.slotStatus) > 0 && iter <= maxIter) || enterTheLoop
                    enterTheLoop = false;
                    % decoding
                    [decRaf,decAcked] = decoding(raf,capturePar);
                    % update acked bursts list
                    % acked.slot   = [acked.slot,decAcked.slot];
                    acked.source = [acked.source,decAcked.source];
                    % perform interference cancellation
                    icRaf = ic(decRaf,sicPar,decAcked);
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

acked.source
queues.status([acked.source])

                % update the packets' status
                output.matrix([sub2ind([input.sources input.queueLength] , transpose(acked.source) , queues.status([acked.source]))]) = output.matrix([sub2ind([input.sources input.queueLength] , transpose(acked.source) , queues.status([acked.source]))]) + 1;

                % TODO: inserire tentativi di ritrasmissione se il pacchetto non è confermato [Issue: https://github.com/afcuttin/jsac/issues/1]
                % update the transmission queues
                % queues.status([acked.source]) = queues.status([acked.source]) - 1;
                % memoryless process (no retransmission attempts)
                queues.status = queues.status - 1;
                source.status = source.status - 1; % update sources statuses
                source.status(source.status < 0) = 0; % idle sources stay idle (see permitted statuses above)

            end

            case 'sul' % random access method is CRDSA

                % TODO: completare Satellite Uplink [Issue: https://github.com/afcuttin/jsac/issues/3]

                load('Captures_SUL','C_SUL');
                capturePar.probability = C_SUL;
                % create the RAF
                numberOfBursts = 2;
                for eachSource1 = 1:source.number
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
                            % decidere se fare qui lo scarto del pacchetto, o più avanti, quando si fa la verifica di quelli confermati
                            % marcare a 0 la posizione corrispondente al pacchetto scartato
                            % decrementare il contatore della coda dei pacchetti in ingresso
                            % procedere con una nuova trasmissione
                            source.status(1,eachSource1)      = 1;
                            [pcktTwins,rafRow]                = generateTwins(raf.length,numberOfBursts);
                            raf.status(eachSource1,pcktTwins) = 1;
                            raf.twins(eachSource1,:)          = rafRow;
                        end
                    end
                end

                acked.slot   = [];
                acked.source = [];

            case {'sdl','tdl'} % no random access, just capture threshold

                % TODO: completare Satellite Downlink [Issue: https://github.com/afcuttin/jsac/issues/6]
                % TODO: completare Terrestrial Downlink [Issue: https://github.com/afcuttin/jsac/issues/5]

                if input.linkMode == 'sdl'
                    load('Captures_SDL','C_SDL');
                    capturePar.probability = C_SDL;
                elseif input.linkMode == 'tdl'
                    load('Captures_TDL','C_TDL');
                    capturePar.probability = C_TDL;
                end

            otherwise
                error('Please select one of the availables link modes (tul, sul, sdl, tdl).');
        end

output.matrix = fliplr(output.matrix); % because packets are updated in reverse order

end