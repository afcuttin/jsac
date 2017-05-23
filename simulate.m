clearvars;

numberOfLoadPoints = 20;
logSpacing         = false;
rafLengths         = [7 10 15];
numberOfSources    = 10;
queueLength        = 1000; % same queue length for every source
% queueLength        = randi([500 700],numberOfSources,1); % variable queue length
queueLength        = 500 % TEST: delete this line after testing
activeModes        = {'tul','sul'}; % a cell array with any of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL

if logSpacing == 1
	proba = logspace(-1,0,numberOfLoadPoints);
elseif logSpacing == 0
	proba = linspace(0.1,1,numberOfLoadPoints);
end

ii = 0;
for modeIndex = activeModes
	for rafIndex = rafLengths
		ii = ii + 1;
		inputPar.rafLength = rafIndex
		for lpIndex = 1:numberOfLoadPoints
			inputPar.poissonThreshold = proba(lpIndex)
			[~,~,~,~,~,~,output] = randomAccess(numberOfSources,queueLength,char(modeIndex),inputPar);
			if numel(queueLength) == 1
				results(ii).load(lpIndex) = queueLength * numberOfSources ./ (output.rafLength * output.duration);
			elseif iscolumn(queueLength)
				results(ii).load(lpIndex) = sum(queueLength) ./ (output.rafLength * output.duration);
			end
			results(ii).throughput(lpIndex) = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
			results(ii).meanRetries(lpIndex) = mean(nonzeros(output.retries));
			assert(results(ii).load(lpIndex) >= results(ii).throughput(lpIndex));
			fprintf('Mode %s - Load %f Throughput %f \n',char(modeIndex),results(ii).load(lpIndex),results(ii).throughput(lpIndex));
			validateResults(queueLength,output);
		end
		results(ii).mode      = char(modeIndex);
		results(ii).rafLength = rafIndex;
	end
end

[path,name,~] = fileparts(mfilename('fullpath'))
datetime = datestr(now,30);
name = strcat('jsac-sim-',datetime);
save(fullfile(path,name),'results');