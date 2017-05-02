# Matlab function to simulate a multiple random access scenario
The matlab function allows to simulate a multiple random access scenario

# Usage
For the purpose of the simulation, use the function as follows:
```
[outQueues,outDelays,outRetries,outFirstTx,outDuration,outRafLength,~] = randomAccess(numberOfSources,queueLength,linkMode);
```
Note the final `~` among the outputs, since that variable (which is a struct) shall be ignored.

# Selection of the modulation scheme and FEC rate
The modulation scheme and the associated FEC code rate are selected using the following:
* [ETSI EN 302 307 V1.3.1](http://www.etsi.org/deliver/etsi_en/302300_302399/302307/01.03.01_20/en_302307v010301a.pdf)
* [DVB-S Bitrate and Bandwidth Calculator](http://www.satbroadcasts.com/DVB-S_Bitrate_and_Bandwidth_Calculator.html)

# Change record
For a list of recent changes, see [here](CHANGES.md)

