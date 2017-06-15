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
	set(groot,'defaultAxesColorOrder',[ 0 0.4470 0.7410; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560 ; 0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330 ; 0.6350 0.0780 0.1840],'defaultAxesLineStyleOrder','-|--|:|-.')

	% plot throughput against load together with a scatterplot of the mean delay
	for figind = unique({results.mode})
		figure('Name',strcat(char(figind),' - Load, Throughput and Delay'),'NumberTitle','off');
		hold on
		title('Load, throughput and delay for different RAF lengths','FontWeight','normal')
		legendInfo = cell(numel(find(strcmp(figind,{results.mode}))),1);
		legendIndex = 1;
		for ploind = find(strcmp(figind,{results.mode}))
			plot(results(ploind).load, results(ploind).throughput);
			legendInfo{legendIndex} = num2str(results(ploind).rafLength);
			legendIndex = legendIndex + 1;
		end
		xlabel('Load')
		ylabel('Throughput')
		colormap jet;
		clrbr = colorbar;
		ylabel(clrbr, 'Mean delay')
		legend(legendInfo,'Location','NorthWest')

		% scatterplot of the mean delay
		for ploind = find(strcmp(figind,{results.mode}))
			scatter(results(ploind).load,results(ploind).throughput,40,(results(ploind).meanDelay),'filled');
		end
	end

	% plot delay against throughput together with a scatterplot of the load
	for figind = unique({results.mode})
		figure('Name',strcat(char(figind),' - Throughput and Delay'),'NumberTitle','off');
		hold on
		title('Throughput, mean delay and load for different RAF lengths','FontWeight','normal')
		legendInfo = cell(numel(find(strcmp(figind,{results.mode}))),1);
		legendIndex = 1;
		for ploind = find(strcmp(figind,{results.mode}))
			plot(results(ploind).throughput,results(ploind).meanDelay);
			legendInfo{legendIndex} = num2str(results(ploind).rafLength);
			legendIndex = legendIndex + 1;
		end
		xlabel('Throughput')
		ylabel('Mean delay')
		clrbr = colorbar;
		ylabel(clrbr, 'Load')
		legend(legendInfo,'Location','SouthEast')

		% scatterplot of the mean load
		for ploind = find(strcmp(figind,{results.mode}))
			scatter(results(ploind).throughput,results(ploind).meanDelay,40,(results(ploind).load),'filled');
		end
	end

	set(groot,'defaultAxesLineStyleOrder','remove')
	set(groot,'defaultAxesColorOrder','remove')
else
	fprintf('Exiting.\n')
end