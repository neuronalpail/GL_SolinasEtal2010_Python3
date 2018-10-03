# GL_SolinasEtal2010_Python3
Cerebellar granular layer model originally built by Solinas et al. 2010^[1]. 

Original model is available at https://github.com/OpenSourceBrain/GranularLayerSolinasNieusDAngelo2010.

The aim of repository is to modify the model for python 3 and to fix bugs.

# Estimated Usage
Type following command in Linux console.

```
$ python GenerateSolinas2010.py
```
	GenerateSolinas2010.py executes nrnivmodl to create x86_64 directory and Start_test.hoc to create SimData directory.

```
$ nrngui Start.hoc
```
	This command will be aborted with error. Start.hoc seems like main simulation file, but there are CVode error.
