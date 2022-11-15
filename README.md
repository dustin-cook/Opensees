# Opensees
This repo can be used to build an Opensees model of a building, run the building to calculate structural response, and post-process the building damage according to the ASCE 41-17 to assess building damage. Scripts are also set up to run the same model through an incremental dynamic analysis (IDA), according to FEMA P-695.

## Important Notes
While the repo is set up to be able to run any building through Opensees, ASCE 41, and IDA, there are currently extensive spots where functionality is hard coded to the Imperial County Services Center. Scripts will need some cleaning up to become general.

## Building TCL Scripts
Build TCL scripts (w/o running any analysis) for all models listed in the models.csv file

NOTE: Currently only set up to build pushover models (different recorders and lateral loads compared to dynamic models)

### Process
- Open driver_build_model_pushover.m
- Edit user inputs section as needed
- Execute script
- Models will populate in the outputs directory


