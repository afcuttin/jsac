% for test purposes only
clearvars;

input.sources     = 100;
input.queueLength = 100;
input.linkMode    = 'tul'; %can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
input.rafLength   = 100;
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)

output = randomAccess(input);

outputLoad       = input.queueLength * input.sources ./ (input.rafLength * output.duration);
outputThroughput = sum(sum(output.matrix,2)) ./ (input.rafLength * output.duration);

fprintf('Results are stored in the "output" variable.\n');
fprintf('Load %f Throughput %f \n',outputLoad,outputThroughput);