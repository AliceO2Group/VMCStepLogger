# MCStepLogger and TMCReplay

### Build and install **with** `aliBuild`

`aliBuild` can be used to setup the package. To obtain `aliBuild` and related packages, please follow https://alice-doc.github.io/alice-analysis-tutorial/building/. To build from the `master` branch, do

```bash
aliBuild init MCStepLogger@master --defaults o2
aliBuild build MCStepLogger --defaults o2
```
or use one of the available [release tags](https://github.com/AliceO2Group/VMCStepLogger/releases).

After that, enter the environment with
```bash
alienv enter MCStepLogger/latest
```
and you are good to go.

## Build and install **without** `aliBuild`

Dependencies are `ROOT` and `Boost` as well as `CMake` in order to build.

`CMake` (version >= 3.15.0) is used to build this project. Please install [ROOT](https://github.com/root-project/root) and [Boost](https://www.boost.org/). For both the versions that should allow the built can be derived from [alidist](https://github.com/alisw/alidist). After that just set the environment variables `ROOTSYS` and `BOOST_ROOT` to the root directories of the `ROOT` and `Boost` installations, respectively. Finally, run
```bash
mkdir -p $BUILD_DIR $INSTALL_DIR; cd $BUILD_DIR
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR $MCSTEPLOGGER_SOURCE_DIR
make install
```
That leaves you with
* libraries at `$INSTALL_DIR/lib`
* headers at`$INSTALL_DIR/include`
* the executable `$INSTALL_DIR/bin/mcStepAnalysis` (usage explained [below](#mcsteploganalysis))

**Additional options**

* In order to disable the built of the `TMCReplay` engine, pass `-MCStepLogger_BUILD_TMCReplay=OFF` at configure time.
* To enable testing, pass `-MCStepLogger_BUILD_TESTS=ON` at configure time.
