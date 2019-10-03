# Voroffset

Discrete mesh offsetting based on half-space Voronoi diagrams and Power diagrams.

### Compilation

Dependencies are downloaded automatically by CMake at configuration time. To compile the code, just type:

```bash
mkdir build
cd build
cmake -j8
```

### Running the code

See possible options with:

```bash
./offset3d -h
```

Example usage:

```
./offset3d filigree.ply -n 512 -p 10 -r 5 -x dilation
```

##### offset2d

Takes a .svg file as input, and saves the result of the dilation as a quad-mesh (.obj).

##### offset3d

Takes a triangle mesh as input (.stl, .obj, .off), and saves the dilated output as a mesh (quad-mesh with .obj, hex-mesh with .mesh, etc.)

```
Offset3D
Usage: ./offset3d [OPTIONS] input [output]

Positionals:
  input TEXT                  Input model
  output TEXT=output.obj      Output model

Options:
  -h,--help                   Print this help message and exit
  -i,--input TEXT             Input model
  -o,--output TEXT=output.obj Output model
  -j,--json TEXT              Output json file
  -d,--dexels_size FLOAT=1    Size of a dexel (in mm)
  -n,--num_dexels INT=256     Number of dexels (-1 to use dexel size instead)
  -p,--padding INT            Padding (in #dexels)
  -t,--num_thread UINT=6      Number of threads
  -r,--radius FLOAT=8         Dilation/erosion radius (in #dexels)
  -m,--method TEXT in {brute_force,ours}
                              The method to use
  -x,--apply TEXT in {closing,dilation,erosion,noop,opening}=dilation
                              Morphological operation to apply
  -f,--force                  Overwrite output file
  -u,--radius_in_mm           Radius is given in mm instead
```

### Replicability

Head over to the [scripts/](scripts/) folder for further instructions.
