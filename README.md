# LinearPaulTrap
My MSc thesis on optimizing a linear Paul trap for single-ion spectroscopy, w/o a buffer gas. The aim is to inject realistic ISOLDE beams inside the trap's apperture, trap them using electric fields and finally cool the ions close to the Doppler limit, using Doppler laser cooling.

The repository includes .gem files for use in SimIon 8.1, User Programs for optimizing the electrodes and analysis codes in Python.

DETAILS:

-> design.gem: The main .gem file to be used by SimIon to create .#PA potential arrays.

-> RFQ_Endcap.lua: The most basic User Program. It relies on applying DC voltages to the Endcap electrodes (No. 10) to slow down and eventually capture the ions.

-> RFQ_Step.lua: Slightly advanced User Program, used to slow down the ions using a step potential.

-> RFQ_Doppler_Cooling.lua: User Program used for Doppler laser cooling, after injecting an ISOLDE beam inside the trap. Some familiarity with the rest of the User Programs is adviced before this one.

-> RFQ_ISOLDE_beam.fly2: User Program that initialises a beam with ISOLDE properties.

-> Analysis_Basic.ipynb and Analysis_Basic_Special.ipynb: Analysis codes to use on raw .txt data obtained from SimIon. The raw data are obtained by printing data according to RFQ_Basic.lua with format: IonNumber & ToF & Energy.

-> Analysis_Step.ipynb: Analysis codes to use on raw .txt data obtained from SimIon. The raw data are obtained by printing data according to RFQ_Step.lua with format: IonNumber & ToF & Energy.

-> Analysis_Doppler_Cooling: Analysis codes to use on raw .txt data obtained from SimIon. The raw data are obtained by:
  - Printing data according to RFQ_Step.lua with format: IonNumber & ToF & Energy
  - Writing data on specific .txt files mentioned in RFQ_Doppler_Cooling.lua
