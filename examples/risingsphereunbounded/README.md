# Rising Sphere Unbounded

## Building
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

## Running
```bash
cd build 
./risingsphereunbounded # Serial
mpirun -np 4 ./risingsphereunbounded
```