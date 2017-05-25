% for test purposes only
clearvars;

numberOfSources            = 10;
queueLength                = 1000; % same queue length for every source
queueLength                = randi([500 700],numberOfSources,1); % variable queue length
inputPars = [];
% inputPars.poissonThreshold = 0.6;

activeModes = {'tul','sul','sdl','tdl'}; % a cell array with any of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL

fprintf('Results are stored in the "results" variable.\n');

for modeIndex = 1:numel(activeModes)
		[~,~,~,~,~,~,output] = randomAccess(numberOfSources,queueLength,char(activeModes(modeIndex)),inputPars);
		if numel(queueLength) == 1
			results(modeIndex).load = queueLength * numberOfSources ./ (output.rafLength * output.duration);
		elseif iscolumn(queueLength)
			results(modeIndex).load = sum(queueLength) ./ (output.rafLength * output.duration);
		end
		results(modeIndex).throughput = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
		results(modeIndex).mode       = char(activeModes(modeIndex));
		fprintf('Mode %s - Load %f Throughput %f \n',results(modeIndex).mode,results(modeIndex).load,results(modeIndex).throughput);
		validateResults(queueLength,output);
end