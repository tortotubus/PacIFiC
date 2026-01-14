# Rising Sphere Unbounded

## Building
```bash
cmake -S . -B build
```

## Running
```bash
cd build 
./risingsphereunbounded # Serial
mpirun -np 4 ./risingsphereunbounded
```