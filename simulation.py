import numpy as np
import abc
from system import *


class Simulation:
    def __init__(self, system, potentials, integrator, outfile='out.pdb'):
        self.system = system
        self.potentials = potentials

        self.time_controller = TimeController()
        self.pdbwriter = PDBWriter(outfile)

        self.euler = Euler(self.potentials, self.time_controller)
        self.integrator = integrator(self.potentials, self.time_controller)

    def start(self, nsteps=100, step=0.01, framerate=20):
        self.pdbwriter.drop_current_state(self.system)
        for i in range(nsteps):
            self.euler.step(self.system, step)
            if i % framerate == 0:
                self.pdbwriter.drop_current_state(self.system)

    def run(self, nsteps=900, step=0.01, framerate=20):
        for i in range(nsteps):
            self.integrator.step(self.system, step)
            if i % framerate == 0:
                self.pdbwriter.drop_current_state(self.system)


class Atom:
    COUNT = 0

    def __init__(self, coords, vel, mass, name):
        self.coords = np.array(coords, dtype='float64')
        self.vel = np.array(vel, dtype='float64')
        self.acc = np.array([0, 0, 0])
        self.mass = mass
        self.name = name
        self.num = Atom.COUNT
        Atom.COUNT += 1

    def move(self, coords):
        self.coords = np.array(coords) #argumentem jest nowa pozycja atomu


class Integrator(abc.ABC):
    @abc.abstractmethod
    def calc_position(self, atom, step):
        pass
    @abc.abstractmethod
    def step(self, system, step):
        pass


class VelocityVerlet(Integrator):

    def __init__(self, potentials, time_controller):
        self.potentials = potentials
        self.time_controller = time_controller

    def calc_position(self, atom, step):
        force = get_force(self.potentials, atom)
        return atom.coords + atom.vel*step + (force*step**2)/(2*atom.mass)

    def step(self, system, step):
        prev = {atom.num: get_force(self.potentials, atom) for atom in system.atoms}
        for atom in system.atoms:
            self.time_controller.add_event(step, atom.move, args=[self.calc_position(atom, step)])
        self.time_controller.exec_step(step)
        for atom in system.atoms:
            force = get_force(self.potentials, atom)
            prevVal = prev[atom.num]
            atom.vel += step*(prevVal+force)/(2*atom.mass)


class Euler(Integrator):

    def __init__(self, potentials, time_controller):
        self.potentials = potentials
        self.time_controller = time_controller

    def calc_position(self, atom, step):
        force = get_force(self.potentials, atom)
        return atom.coords + atom.vel*step + (force*step**2)/(2*atom.mass)

    def step(self, system, step):
        for atom in system.atoms:
            self.time_controller.add_event(step, atom.move, args=[self.calc_position(atom, step)])
        self.time_controller.exec_step(step)
        for atom in system.atoms:
            force = get_force(self.potentials, atom)
            atom.vel += force*step/atom.mass


def get_force(potentials, atom):
    force = np.zeros(3)
    for potential in potentials:
        force += potential.calc_force(atom)
    return np.array(force, dtype='float64')
