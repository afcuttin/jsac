# Matlab function to simulate a multiple random access scenario
The matlab function allows to simulate a multiple random access scenario

# Change record
For a list of recent changes, see [here](CHANGES.md)

# Usage
For the purpose of the simulation, use the function as follows:
```
[outQueues,outDelays,outRetries,outFirstTx,outDuration,outRafLength,~] = randomAccess(numberOfSources,queueLength,linkMode);
```
Note the final `~` among the outputs, since that variable (which is a struct) shall be ignored.

# Simulation parameters
The following are the parameter to be used for the joint simulation.

In every scenario:
* the maximum number of available iterations for the decoding and the interference cancellation is 1;
* the maximum number of possible retransmission (retry limit) in case a packet is not acknowledged is 4.
* the maximum number of active sources must not be grater than 10.

The modulation parameters are as follows:
* 8-PSK modulation
* the FEC code rate is 2/3 or 1/3

<!---
## Terrestrial uplink (tul) scenario

## Satellite uplink (sul) scenario

## Satellite downlink (sdl) scenario

## Terrestrial downlink (tdl) scenario
-->

# Selection of the modulation scheme and FEC rate
The modulation scheme and the associated FEC code rate are selected using the following:
* [ETSI EN 302 307 V1.3.1](http://www.etsi.org/deliver/etsi_en/302300_302399/302307/01.03.01_20/en_302307v010301a.pdf)
* [DVB-S Bitrate and Bandwidth Calculator](http://www.satbroadcasts.com/DVB-S_Bitrate_and_Bandwidth_Calculator.html)

# Usage of .mat files
* N_v: vettore del numero di sorgenti totali (desiderato +interferenti) nello slot (da 1 a 10)
* S_v: vettore delle soglie (in dB) sul SINR (da -10 a 10)
* R_v: vettore delle soglie sul tasso - nel senso di efficienza spetrale (modulazione e tasso del codice) (da 0.05 a 10)

## Captures TUL
* C\_TUL\_3(numero di sorgenti nel primo slot,numero di sorgenti nel secondo slot,numero di sorgenti nel terzo slot,indice soglia sul tasso)
* C\_TUL\_4(numero di sorgenti nel primo slot,numero di sorgenti nel secondo slot,numero di sorgenti nel terzo slot,numero di sorgenti nel quarto slot,indice soglia sul tasso)

## Captures SUL
* C_SUL(indice numero di sorgenti,indice soglia sul SINR)

## Captures SDL
* C_SDL(indice soglia sul SINR)

## Captures TDL
* C_TDL(indice soglia sul SINR)