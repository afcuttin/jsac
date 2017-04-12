% for test purposes only
clearvars;

numberOfSources           = 20;
queueLength               = 100;
% queueLength               = transpose(randperm(numberOfSources*2,numberOfSources));

activeModes = {'sul'}; % a cell array with any of the following: 'tul', terrestrial uplink, 'sul', satellite UL, 'sdl', satellite downlink, 'tdl', terrestrial DL

fprintf('Results are stored in the "output" variable.\n');

% TODO: every link mode shall be tested [Issue: https://github.com/afcuttin/jsac/issues/33]
if any(strcmp('tul',activeModes))

	output         = randomAccess(numberOfSources,queueLength,'tul');
	% TODO: il calcolo del throughput deve essere differenziato a seconda che si usi CRDSA (slot) o CSA (slice). Per adesso sono uguali [Issue: https://github.com/afcuttin/jsac/issues/12]
	TUL_Load       = queueLength * numberOfSources ./ (output.rafLength * output.duration);
	TUL_Throughput = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
	fprintf('TUL - Load %f Throughput %f \n',outputLoad,outputThroughput);

elseif any(strcmp('sul',activeModes))

	output         = randomAccess(numberOfSources,queueLength,'sul');
	SUL_Load       = queueLength * numberOfSources ./ (output.rafLength * output.duration);
	SUL_Throughput = sum(sum(output.queues,2)) ./ (output.rafLength * output.duration);
	fprintf('SUL - Load %f Throughput %f \n',outputLoad,outputThroughput);

elseif any(strcmp('tul',activeModes))

	output = randomAccess(numberOfSources,queueLength,'sdl');
	% SUL_Load
	% SUL_Throughput
	fprintf('SDL - Load %f Throughput %f \n',outputLoad,outputThroughput);

elseif any(strcmp('tul',activeModes))

	output = randomAccess(numberOfSources,queueLength,'tdl');
	% SUL_Load
	% SUL_Throughput
	fprintf('TDL - Load %f Throughput %f \n',outputLoad,outputThroughput);

end
