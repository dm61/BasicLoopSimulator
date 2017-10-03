# BasicLoopSimulator 
Script runs in MATLAB (https://www.mathworks.com/) or in Octave (https://www.gnu.org/software/octave/)

Assumptions and limitations 
* Simulation of a single meal with known carb absorption  
* Initially IOB(0)=0, COB(0)=0
* Ideal system model: deltaBG = (carb impact)-(insulin impact)
* Does not incude Loop retrospective correction or momentum effects, or dynamic carb algorithm
* Assumes a single fixed nominal basal rate during entire simulation
* Includes exponential insulin absorption curves with td, tp parameters
* Includes Loop v1.4 and Loop v1.5 dosing algorithms
