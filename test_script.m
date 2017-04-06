% for test purposes only
clearvars;

input.sources             = 200;
input.queueLength         = 300;
input.linkMode            = 'tul'; %can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
input.rafLength           = 80;
input.burstMaxRepetitions = 4;
input.bitsPerSymbol       = 2;
input.fecRate             = 4/5;
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)

output = randomAccess(input);

if strcmp(input.linkMode,'tul')
	% TODO: il calcolo del throughput deve essere differenziato a seconda che si usi CRDSA (slot) o CSA (slice). Per adesso sono uguali [Issue: https://github.com/afcuttin/jsac/issues/12]
	outputLoad       = input.queueLength * input.sources ./ (input.rafLength * output.duration);
	outputThroughput = sum(sum(output.matrix,2)) ./ (input.rafLength * output.duration);
elseif strcmp(input.linkMode,'tdl') || strcmp(input.linkMode,'sdl') || strcmp(input.linkMode,'sul')
	outputLoad       = input.queueLength * input.sources ./ (input.rafLength * output.duration);
	outputThroughput = sum(sum(output.matrix,2)) ./ (input.rafLength * output.duration);
end

fprintf('Results are stored in the "output" variable.\n');
fprintf('Load %f Throughput %f \n',outputLoad,outputThroughput);
