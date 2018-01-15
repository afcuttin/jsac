% for test purposes only
clearvars;

numberOfSources = 5;
queueLength     = randi([500 700],numberOfSources,1); % variable queue length
queueLength     = 5000; % same queue length for every source
inputPars       = [];
activeModes     = {'tul','sul','sdl','tdl'}; % a cell array with any of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL

% inputPars.retryLimit       = 4;
% inputPars.activeModes      = {'tul'};
% inputPars.sicMaxIter       = 2;
% inputPars.tulAccMeth       = 'csa-nc';
% inputPars.tulAccMeth       = 'csa';
% inputPars.poissonThreshold = 0.0116;
% inputPars.poissonThreshold = 0.75;
% inputPars.rafLength = 500;
% inputPars.rafLength = 10;
% inputPars.timeConstraint  = 1000;
% inputPars.timeConstraint  = 10;

fprintf('Results are stored in the "results" variable.\n');

for modeIndex = 1:numel(activeModes)
		[~,~,~,~,~,~,output]          = randomAccess(numberOfSources,queueLength,char(activeModes(modeIndex)),inputPars);
		results(modeIndex).load       = sum(sum(output.firstTx > 0,2)) ./ (output.rafLength * output.duration);
		results(modeIndex).throughput = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
		results(modeIndex).mode       = char(activeModes(modeIndex));
		fprintf('Mode %s - Load %f Throughput %f \n',results(modeIndex).mode,results(modeIndex).load,results(modeIndex).throughput);
		validateResults(queueLength,output);
end