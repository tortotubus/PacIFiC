# Apptainer

Build any .def files into portable images using 
```bash
apptainer build mycontainer.sif mycontainer.def
```
these containers then may be used interactively with `apptainer shell mycontainer.sif` or have commands run inside them using `apptainer exec mycontainer.sif my_command`. 

The hybrid container definitions build the defined version of MPI and will work with the host cluster's `mpirun` so long as the MPI versions in both the container and on the cluster are ABI-compatible.

These images may be built on any machine and then copied onto the cluster with the caveats stated above. It is advisable to use the hybrid containers if possible since cluster installations of MPI are often optimized for their own hardware/interconnects.

Apptainer Docs:
* https://apptainer.org/docs/user/latest/mpi.html

HPC Cluster Docs:
* https://docs.alliancecan.ca/wiki/Apptainer
* https://confluence.it.ubc.ca/spaces/UARC/pages/279849945/Using+Apptainer+or+Singularity+Containers

