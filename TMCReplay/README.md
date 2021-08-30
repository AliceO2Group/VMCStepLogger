# TMCReplay

This is a pseudo-detector-simulation engine based on the [Virtual Monte Carlo (VMC)](https://vmc-project.github.io/) package. It takes logged steps recorded by the `MCStepLogger` to **replay** the entire simulation.

The main objective is to be able to provide an engine to study the impact of parameter variations with. Eventually, this can be used to optimise (full-)simulation parameters in view of enhancing their speed and efficiency.

The functionality is compiled into a separate library `libTMCReplay` which also allows to use the `libMCStepLoggerIntercept` which also allows to run the step logging against the `TMCReplay` engine.

## Setting parameters

The list of tunable parameters is foreseen to correspond to the possible settings implemented for [GEANT3_VMC]() and [GEANT4_VMC](https://github.com/vmc-project/geant4_vmc) and can be found [here](include/TMCReplay/Physics.h). The parameters can be set via the usual interfaces
* `TVirtualMC::SetCut` and `TVirtualMC::SetProcess` to set global cut and process values
* `TVirtualMC::Gstpar` to set parameters depending on a certain medium

**Please note** that the parameters can be set but at the moment nothing is actually applied. Therefore, the exact same steps would be reproduced regardless of the parameters. Of course, this is the next major development as it is one of the core objectives of this engine.

## Testing different parameter settings

**Note that this section becomes actually meaningful when the aforementioned functionality will have been implemented.** However, in the current state it can already be verified that the `TMCReplay` engine does what it is supposed to do.

### Running/testing on step level only (`TMCReplayDummyApplication` and `TMCReplayDummyStack`)

If only the impact on the steps is desired, all required functionality is already included here. Simply have a look at [this executable](src/replay.cxx). It makes use of the `TMCReplayDummyApplication` and `TMCReplayDummyStack` just to replay the steps. Once implemented, the impact of the parameters on the stepping can be evaluated with this tool alone. To do so, just run the `MCStepLogger` together with the replay, for instance as follows

```bash
MCSTEPLOG_TTREE=1 LD_PRELOAD=$MCSTEPLOGGER_ROOT/lib/libMCStepLoggerIntercept.so tmcreplay
```

That will pick up the step file `MCStepLoggerOutput.root`. If you gave that a different name (or in case also your contained tree has another name), you will find the following help message useful:

```bash
> tmcreplay --help
Simulation options:
  --help                                show this help message and exit
  --stepfilename arg (=MCStepLoggerOutput.root)
                                        MCStepLogger filename
  --steptreename arg (=StepLoggerTree)  treename inside file where to find step
                                        tree
  --geofilename arg (=o2sim_geometry.root)
                                        ROOT geometry filename
  --geokeyname arg (=FAIRGeom)          key name inside geo file where to find
                                        geometry tree
```

As you can see the geometry is required by this executable as well since it has otherwise no idea what geometry was used for the original simulation.

### Running in another environment

In that case the geometry construction is assumed to implemented by that environment's implementation of `TVirtualMCApplication`. Also, no stack is constructed automatically and it has to be constructed and passed by the user.

Such a scenario is useful to evaluate the impact of parameter variation on hits or even digits.
