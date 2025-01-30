To compile the package:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc11-opt/setup.sh
mkdir build
cd build
cmake ..
make
```

To run the package (from the `tools` folder):
```
build/apps/angles.exe Lb2L0mm_rapidsim.json
build/apps/angles.exe Lb2L0mm_reconstructed.json
```
You will need to change the output location to someplace where you have write access.
Also, don't trust the constrained variables, this bit does not work.