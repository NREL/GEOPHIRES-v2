# GEOPHIRES v2.0
## GEOPHIRES v2.0 Description
GEOPHIRES v2.0 is a free and open-source geothermal techno-economic simulator. GEOPHIRES combines reservoir, wellbore, surface plant, and economic models to estimate the capital and operation and maintenance costs, instantaneous and lifetime energy production, and overall levelized cost of energy of a geothermal plant. Various reservoir conditions (EGS, doublets, etc.) and end-use options (electricity, direct-use heat, cogeneration) can be modeled. The tool can be coupled to the external reservoir simulation TOUGH2. Users are encouraged to build upon to the GEOPHIRES v2.0 framework to implement their own correlations and models.

## GitHub Folder Info
This GitHub folder contains the following folders and files:
- `GEOPHIRESv2.py`: GEOPHIRES v2 Python code
- `GEOPHIRES v2.0 User Manual.pdf`: User manual including quick start guide and list of all input parameters.
- `References/`: Folder containing reference documents on GEOPHIRES
- `Examples/`: Folder containing example problems

## Importing as Python Module

Example:

```python
import os
GEOPHIRES_v2 = __import__('GEOPHIRES-v2.GEOPHIRESv2')

if __name__ == '__main__':
    GEOPHIRES_v2.GEOPHIRESv2.run_geophires(input_file_path=os.path.join('data', 'my-geophire-inputs.txt'))
```

## Contact
In case of questions, comments, or suggestions for improvement or collaboration, please contact Koenraad Beckers (koenraad.beckers@heateon.com) or Kevin McCabe (kevin.mccabe@nrel.gov).

## License
GEOPHIRES v2.0 is distributed as open-source under the GNU General Public License version 3.0 (GPLv3), Copyright (c) 2019 Alliance for Sustainable Energy, LLC.

## Version
- GEOPHIRES v2.0 Build Date: December 9, 2018