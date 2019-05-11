# EventWeightReader

This project is what is used to generate the NuMI nue CC inclusive cross section systematics. The main function is the larsoft module named `EventWeightReader_module.cc`. This module will take an artroot file which contains the mc weight data product, get the weights and return a data structure containing the reweighted events in each universe (whether it will be a unisim or multisim). There are various ingredients for this needed in order for this to work. The first thing is you will need a list of artroot files from the output of the selection that have been ran through `EventWeight`. For the cross section re-calculation, these need to have both the generated in tpc signal events and any selected events. You do not need the background events that dont make the selection. 

You also need a txt file which contains the event numbers of the reweighted events. This is so the module can reweight the correct events. (See the module code to see what format the txt files need to be in, its fairly felxible depending on the users need). 

The output of the module is a root file which contains an Event_list class which contains all the information needed to get the systematic uncertainty. Once you have this file, run the make_xsection_plot.C root macro to take a look at the systematics. 
  
Before running, you will need to create the event_list library by running make from the `make` file and building larsoft as usual  

commands to run
`lar -c run_EventWeightReader.fcl /uboone/data/users/kmistry/work/NueXSection_Outputs/eventweight/v2/genie/filtered_eventweight_genie_multisim_250Univ_individual.root`   
`lar -c run_EventWeightReader.fcl /uboone/data/users/kmistry/work/NueXSection_Outputs/eventweight/v2/model/filtered_eventweight_model_multisim_1000Univ_individual.root `   
`lar -c run_EventWeightReader.fcl /uboone/data/users/kmistry/work/NueXSection_Outputs/eventweight/v2/reinteractions/filtered_eventweight_reinteractions_multisim_250Univ_individual.root`  


I have not designed this for other people to use, but if anyone was intersted in using the code, then get in contact with me!
