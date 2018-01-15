clearvars;

% input parameters for the simulation:
numberOfLoadPoints       = 40;
logSpacing               = false;
% rafLengths               = [500];
rafLengths               = [35];
% inputPar.numberOfSources = 20000;
inputPar.numberOfSources = 100;
inputPar.queueLength     = 1000000; % same queue length for every source
% inputPar.queueLength   = randi([1800 2200],numberOfSources,1); % variable queue length
% inputPar.activeModes     = {'sul'}; % a cell array with any of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL
inputPar.activeModes     = {'tul'}; % a cell array with any of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL
inputPar.timeConstraint  = 100000;
inputPar.retryLimit      = 1;
inputPar.sicMinIter      = 1;
% inputPar.sicMaxIter      = inf(1);
inputPar.sicMaxIter      = 2;
% inputPar.tulAccMeth      = 'csa-nc';
inputPar.tulAccMeth      = 'csa';
% inputPar.sulAccMeth      = 'crdsa-nc';

% Goff = logspace(-2,0,numberOfLoadPoints)
Goff = linspace(0.01,1.5,numberOfLoadPoints)
% Goff = [0.01:0.05:0.4 0.4:0.02:0.6 0.6:0.05:1]
% numberOfLoadPoints = numel(Goff);

if logSpacing == 1
	proba = logspace(-3,-0.7,numberOfLoadPoints);
elseif logSpacing == 0
	proba = linspace(0.01,0.017,numberOfLoadPoints); %csa
	proba = linspace(0.001,0.025,numberOfLoadPoints); % crdsa
	proba = Goff * rafLengths / inputPar.numberOfSources
	proba(proba>1) = 1;
end

fprintf('Started at \n')
fprintf('%u ',fix(clock))
fprintf('\n')
t1 = tic;
ii = 0;
for modeIndex = inputPar.activeModes
	for rafIndex = rafLengths
		ii = ii + 1;
		inputPar.rafLength = rafIndex;
		for lpIndex = 1:numberOfLoadPoints
            t2=tic;
			inputPar.poissonThreshold = proba(lpIndex);
			[~,~,~,~,~,~,output] = randomAccess(inputPar.numberOfSources,inputPar.queueLength,char(modeIndex),inputPar);
			% if numel(queueLength) == 1
				% results(ii).load(lpIndex) = queueLength * numberOfSources ./ (output.rafLength * output.duration);
			% elseif iscolumn(queueLength)
				results(ii).load(lpIndex) = sum(sum(output.firstTx > 0,2)) ./ (output.rafLength * output.duration);
			% end
			results(ii).throughput(lpIndex) = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
			results(ii).meanDelay(lpIndex) = mean(nonzeros( ((output.retries -1) .* output.rafLength + output.delaySlot) .* output.queues ));
			assert(results(ii).load(lpIndex) >= results(ii).throughput(lpIndex));
			validateResults(inputPar.queueLength,output);
			fprintf('Mode %s - Load %f Throughput %f. \n %.0f seconds for Simulation %u of %u. Elapsed time %u m %.0f s. \n',char(modeIndex),results(ii).load(lpIndex),results(ii).throughput(lpIndex),toc(t2),(lpIndex + (ii-1)*numel(proba)),numel(inputPar.activeModes)*numel(rafLengths)*numberOfLoadPoints,floor(toc(t1)/60),toc(t1)-floor(toc(t1)/60)*60);
		end
		results(ii).mode      = char(modeIndex);
		results(ii).rafLength = rafIndex;
	end
end

[path,~,~] = fileparts(mfilename('fullpath'));
datetime = datestr(now,30);
name = strcat('jsac-sim-',datetime);
save(fullfile(path,name));
totalTime = toc(t1);
fprintf('Total time elapsed: %u hours %2.0f minutes (%.0f seconds).\n',fix(totalTime/3600),(fix(totalTime/60)-fix(totalTime/3600)*60),totalTime);
fprintf('Ended at \n')
fprintf('%u ',fix(clock))
fprintf('\n')