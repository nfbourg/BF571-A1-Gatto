# BF571-A1-Gatto
Move this file into the Matlab folder. 
It is important to install the Raventoolbox and Mosek.

The file named "getSGDinRefModels" is the model from the Gatto paper. 
The following commands will reproduce figures from the Gatto paper:
getSGDinRefModels('iRC1410','FBS','woFluxes')
getSGDinRefModels('iRC1410','HAM','woFluxes')
getSGDinRefModels('iPC1675','HAM','woFluxes')

The custom model used for our project is run by first importing the models using the runFBAs script.
Then you will run the projectfba script after editing in the model name you want to use and the setting the respective optimization parameter at the top of the script.
