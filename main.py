from potentials import *
from simulation import *
from system import *
import argparse
from random import randint
from random import sample

parser = argparse.ArgumentParser(description="Make animation of the particles in the box (return PDB file)")
parser.add_argument('--atom', '-a', default=['He', 'C'], metavar='SYMBOL', type=str, action='append',
                    help='add atom to the box')
parser.add_argument('--bonded', '-b', default=[['H', 'H'], ['H', 'H']], metavar='SYMBOL', nargs=2, type=str, action='append',
                    help='add two bonded particles')
parser.add_argument('--atoms', '-aa', default=[[2, 'O']], metavar='NUMBER SYMBOL', nargs=2, action='append',
                    help='add number of particles of the same type to the box')
parser.add_argument('--box', '-x', default=[10, 10, 10], metavar='WIDTH LENGTH HEIGHT', nargs=2, action='store',
                    help='add box of the given size')
parser.add_argument('--steps', '-s', default=10000, metavar='NUMBER', action='store',
                    help='number of steps of the simulation')
parser.add_argument('--output', '-o', type=str, default='out.pdb', metavar='FILE_NAME', action='store',
                    help='name of the output pdb file')

args = parser.parse_args()

weights = {'H': 1, 'He': 4, 'Li': 7, 'Be': 9, 'B': 10, 'C': 12, 'N': 14, 'O': 16, 'F': 19, 'Na': 23, 'Mg': 24, 'Al': 27,
           'Si': 28, 'P': 31, 'S': 32, 'Cl': 35, 'K': 39, 'Ca': 40}
vel = [-1, 1]*2 + [0]*4

system = System()
potentials = [VanDerWaals(system)]
if args.box != [0, 0, 0]:
    w, l, h = args.box
    potentials += [Wall([-1, 0, 0], d=w/2), Wall([1, 0, 0], d=w/2)]
    potentials += [Wall([0, -1, 0], d=l/2), Wall([0, 1, 0], d=l/2)]
    potentials += [Wall([0, 0, -1], d=h/2), Wall([0, 0, 1], d=h/2)]
else:
    w, l, h = 10, 10, 10

for a in args.atom:
    system.add(Atom([randint(-w/2+1, w/2-1), randint(-l/2+1, l/2-1), randint(-h/2+1, h/2-1)], sample(vel, 3), weights[a], a))

for bond in args.bonded:
    location = [randint(-w/2+1, w/2-1), randint(-l/2+1, l/2-1), randint(-h/2+1, h/2-1)]
    velocity = sample(vel, 3)
    a1 = Atom(location, velocity, weights[bond[0]], bond[0])
    a2 = Atom([l+ll for l, ll in zip(location, [1, 0, 0])], velocity, weights[bond[1]], bond[1])
    system.add(a1)
    system.add(a2)
    potentials.append(Bond(a1, a2))

for num, a in args.atoms:
    for _ in range(num):
        system.add(Atom([randint(-w/2+1, w/2-1), randint(-l/2+1, l/2-1), randint(-h/2+1, h/2-1)], sample(vel, 3), weights[a], a))

simulation = Simulation(system, potentials, VelocityVerlet, outfile=args.output)
simulation.start()
simulation.run(nsteps=args.steps*2)
simulation.pdbwriter.close()
