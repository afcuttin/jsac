% for test purposes only
clear all;

input.sources     = 400;
input.queueLength = 1000;
input.linkMode    = 'tul'; %can be one of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL (type: string)
input.rafLength   = 100;
% * input.sinrThreshold: value of the SINR threshold to be used (type: integer)

output = randomAccess(input);

outputLoad       = input.queueLength * input.sources ./ (input.rafLength * output.duration)
outputThroughput = sum(sum(output.matrix,2)) ./ (input.rafLength * output.duration)