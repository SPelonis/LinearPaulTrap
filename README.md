# LinearPaulTrap
My MSc thesis on optimizing a linear Paul trap for single-ion spectroscopy, w/o a buffer gas. The aim is to inject realistic ISOLDE beams inside the trap's apperture, trap them using electric fields and finally cool the ions close to the Doppler limit, using Doppler laser cooling.

The repository includes .gem files for use in SimIon 8.1, User Programs for optimizing the electrodes and analysis codes in Python.

DETAILS:

-> design.gem: The main .gem file to be used by SimIon to create .#PA potential arrays.

-> RFQ_Endcap.lua: The most basic User Program. It relies on applying DC voltages to the Endcap electrodes (No. 10) to slow down and eventually capture the ions.

-> RFQ_Step.lua: Slightly advanced User Program, used to slow down the ions using a step potential.

-> RFQ_Doppler_Cooling.lua: User Program used for Doppler laser cooling, after injecting an ISOLDE beam inside the trap. Some familiarity with the rest of the User Programs is adviced before this one.
