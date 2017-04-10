% for test purposes only
clearvars;

numberOfSources           = 80;
queueLength               = transpose(randperm(numberOfSources*2,numberOfSources));
linkMode                  = 'tul'; %can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)

% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)

output = randomAccess(numberOfSources,queueLength,linkMode);

if strcmp(linkMode,'tul')
	% TODO: il calcolo del throughput deve essere differenziato a seconda che si usi CRDSA (slot) o CSA (slice). Per adesso sono uguali [Issue: https://github.com/afcuttin/jsac/issues/12]
	outputLoad       = queueLength * numberOfSources ./ (output.rafLength * output.duration);
	outputThroughput = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
elseif strcmp(linkMode,'tdl') || strcmp(linkMode,'sdl') || strcmp(linkMode,'sul')
	outputLoad       = queueLength * numberOfSources ./ (output.rafLength * output.duration);
	outputThroughput = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
end

fprintf('Results are stored in the "output" variable.\n');
fprintf('Load %f Throughput %f \n',outputLoad,outputThroughput);
