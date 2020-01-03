# Investigating `Everything`


The following minimal program can be used.

```cpp
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/Util.h>
#include <commands/parser.h>

void debug_Everything( Everything *e )
{
    // print out mysterious variables here ...
}

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);

    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();
    debug_Everything(&e);

    return 0;
}
```

Build the program using `build.sh` script.

To minimize the time required by `e.setup()`, put the following line in the input file.
```
wavefunction random
```

# Wavefunctions

Bloch wavefunctions are represented as `vector<ColumnBundle>` in `Everything.eVars`.
It can be accessed by: `e.eVars.C`.

The length of `e.eVars.C` can be obtained using `.size()` method.
It is the same as `Nkspin` in `PWDFT.jl`. In JDFTX it is named `nStates`.

One wave function in a k-point is represented by `ColumnBundle`.
Number of electronic states can be accessed using `.nCols()` method.
Number of basis function can be accessed using `.colLength()` method.

# Calculating electron density

```c++
//! Calculate density using current orthonormal wavefunctions (C)
ScalarFieldArray calcDensity() const;
```

# Electronic minimization

For my normal use case IonicMinize is used.

It seems that ElecMinimize can not be used directly, at least for metallic case.
For metallic case, the function `elecFluidMininimize()` should be used.

```cpp
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <commands/parser.h>

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);
    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();

    elecFluidMinimize(e);

    return 0;
}
```