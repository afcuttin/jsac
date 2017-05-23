% make the plots in this file
resultsExist = exist('results','var');
proceed = 0;
if resultsExist == 1
	prompt = 'Do you want to load results different from the existing ones? y/n [n]: ';
	str = input(prompt,'s');
	if isempty(str)
	    str = 'n';
	    fprintf('I will plot existing result.\n');
	    proceed = 1;
	elseif str == 'y'
		prompt = 'This will erase existing the results variable. Are you sure? y/n [n]: ';
		str = input(prompt,'s');
		if isempty(str) || str == 'n'
	    	str = 'n';
	    	fprintf('I will plot existing result.\n');
	    	proceed = 1;
		elseif str == 'y'
			clear('results');
			[filename, ~,~] = uigetfile('jsac-sim*.mat', 'Select a .mat file');
			if ~isempty(filename)
				proceed = 1;
				load(filename);
			end
		end
	else
		fprintf('Please answer "y" or "n". Start again.\n');
	end
elseif resultsExist == 0
	fprintf('There are no results to plot. Looking for a file to load them.\n');
	[filename, ~,~] = uigetfile('jsac-sim*.mat', 'Select a .mat file');
	if filename ~= 0
		proceed = 1;
		load(filename);
	end
end

if proceed == 1
	set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],'defaultAxesLineStyleOrder','-|--|:')
	% plot throughput against load together with a scatterplot of the mean delay
	for figind = unique({results.mode})
		figure('Name',char(figind),'NumberTitle','off');
		hold on
		for ploind = find(strcmp(figind,{results.mode}))
			plot(results(ploind).load, results(ploind).throughput);
			scatter(results(ploind).load,results(ploind).throughput,40,(results(ploind).meanRetries),'filled');
		end
	end

	set(groot,'defaultAxesLineStyleOrder','remove')
	set(groot,'defaultAxesColorOrder','remove')
else
	fprintf('Exiting.\n')
end