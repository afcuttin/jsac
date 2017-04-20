function [response] = validateResults(queueLength,input)

for ii=1:numel(queueLength)
	assert(sum((input.queues(ii,[1:1:queueLength(ii)]) == 0),2) + sum(input.queues(ii,:),2) == queueLength(ii),'Test failed (row %u)\n',ii)
end

if all(all(input.queues == (input.delays == input.firstTx .* input.queues + (input.retries - 1)))) == 1
	response = 'passed';
elseif all(all(input.queues == (input.delays == input.firstTx .* input.queues + (input.retries - 1))))  ~= 1
	response = 'failed';
end

fprintf('Test %s \n',response);