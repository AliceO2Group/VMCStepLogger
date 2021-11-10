# MCReplay

This is a pseudo-detector-simulation engine based on the [Virtual Monte Carlo (VMC)](https://vmc-project.github.io/) package. It takes logged steps recorded by the `MCStepLogger` to **replay** the entire simulation.

The main objective is to be able to provide an engine to study the impact of parameter variations with. Eventually, this can be used to optimise (full-)simulation parameters in view of enhancing their speed and efficiency.

The functionality is compiled into a separate library `libMCReplayCore`. Hence, `libMCStepLoggerInterceptSteps` is completely independent which also allows to run the step logging against the `MCReplayEngine` engine.

The replay has been tested and is verified against [GEANT3_VMC](https://github.com/vmc-project/geant3) and [GEANT4_VMC](https://github.com/vmc-project/geant4_vmc) VMC interfaces.

## Setting parameters

The list of tunable parameters is foreseen to correspond to the possible settings implemented for GEANT3_VMC and GEANT4_VMC and can be found [here](include/MCReplay/MCReplayPhysics.h). The parameters can be set via the usual interfaces
* `TVirtualMC::SetCut` and `TVirtualMC::SetProcess` to set global cut and process values,
* `TVirtualMC::Gstpar` to set parameters depending on a certain medium.
The cut parameters act like production cuts.

**Please note** that the parameters can be set but at the moment nothing is actually applied. Therefore, the exact same steps would be reproduced regardless of the parameters. Of course, this is the next major development as it is one of the core objectives of this engine.
**On the other hand**, there is one special parameter called `CUTALLE` which is an energy production cut on **any** PDG which can also be set through the aforementioned interfaces. Note, that this parameter has no meaning to GEANT3_VMC or GEANT4_VMC.

## Testing different parameter settings

**Note that this section becomes actually meaningful when the aforementioned functionality will have been implemented.** However, in the current state it can already be verified that the `MCReplayEngine` engine does what it is supposed to do.

### Running/testing on step level only (`MCReplayGenericApplication` and `MCReplayGenericStack`)

If only the replay of steps is desired, all required functionality is already included here, completely independent of any environment the reference run was done. Simply have a look at [this code](src/replay.cxx) serving as the source for the executable `mcreplay`. It makes use of the `MCReplayGenericApplication` and `MCReplayGenericStack` just to replay the steps. Once implemented, the impact of the parameters on the stepping can be evaluated with this tool alone. To do so, just run the `MCStepLogger` followed by the replay, for instance as follows (assuming  everything was setup with `aliBuild`)

```bash
MCSTEPLOG_TTREE=1 LD_PRELOAD=$MCSTEPLOGGER_ROOT/lib/libMCStepLoggerIntercept.so mcreplay
```

That will pick up the step file `MCStepLoggerOutput.root`. If you gave that a different name (or in case also your contained tree has another name), you will find the following help message useful:

```bash
> mcreplay --help
Replaying a previously recorded particle transport step-by-step. Mainly meant for checking and performance measurements/optimisation of the transport:
  --help                                show this help message and exit
  --stepfilename arg (=MCStepLoggerOutput.root)
                                        MCStepLogger filename
  --steptreename arg (=StepLoggerTree)  treename inside file where to find step
                                        tree
  --geofilename arg (=o2sim_geometry.root)
                                        ROOT geometry filename
  --geokeyname arg                      key name inside geo file where to find
                                        geometry tree
  -n [ --nevents ] arg (=-1)            number of events to replay
  -e [ --energycut ] arg (=-1)          energy cut to be applied [GeV]
```

As you can see the geometry is required by this executable as well since it has otherwise no idea what the geometry should look like from the `MCStepLoggerOutput.root` file alone.

### Running in another environment

In that case the geometry construction is assumed to implemented by that environment's implementation of `TVirtualMCApplication`. Also, no stack is constructed automatically and it has to be constructed and passed by the user.

Such a scenario is useful to evaluate the impact of parameter variation on hits or even digits.

## Using the engine as a full simulation engine

The `MCReplayEngine` has been tested with the [ALICE O2 framework](https://github.com/AliceO2Group/AliceO2) and also hits can be reproduced sufficiently which allows for realistic parameter tuning.
