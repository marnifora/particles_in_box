# Particles in a box

### Simulation of particles' movement in a limited space

<br></br>

#### Versions:
Python: 3.6.6

numpy: 1.16.3

<br></br>


#### How to use?

./python3 main.py

See main.py --help for more information about changing parameters such as number of atoms, bonds between them or size of the box.

Default number of atoms is 2 oxygens, 1 carbon, 1 helium and 2 diatomic hydrogen molecules. 
Default size of the box is 10x10x10, standard number of steps is 10 000.

Output of the script is pdb file with number of states (as film frames). 
It can be easily displayed using for example Pymol.
