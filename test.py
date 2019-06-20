from simulation import *
from potentials import *
from system import *
import numpy as np


def builder1():

    atoms = [Atom(np.array([0.0,0.0,0.0]), 1, 'a1', np.array([0.0,0.0,1.0])),
             Atom(np.array([0.0,0.0,10.0]), 1, 'a2', np.array([0.0,0.0,-1.0]))]
    s = System()
    for a in atoms:
        s.add(a)
    return s


def test1(t=10):

    s = builder1()
    writer = PDBWriter('test1.pdb')
    time = TimeController()
    print('atomy systemu'+str(s.atoms))
    for i in range(t):
        writer.drop_current_state(s)
        for a in s.atoms:
            time.add_event(1, a.move, [a])
        print('aktualny czas: '+str(i))
        time.exec_step(i)


def builder3(n=10):

    s = System()
    for i in range(n):
        coords = np.random.choice(np.arange(20), 3)
        vel = np.random.choice(np.arange(0, 1, 0.1), 3)
        mass = 1
        s.add(Atom(coords, mass, 'a%s' % i, vel))
    return s


def test3(t=10):

    s = builder1()  
    time = TimeController()
    integrator = VelocityVerlet()
    euler = Euler()
    p = np.array([VanDerWaals()])

    s = Simulation(time, s, integrator, euler, p)

    writer = PDBWriter('test4.pdb')
    s.compile(t, writer)


def test_wall(k=4):

    nor = np.array([1,1,0])
    d = 1
    wall1 = Wall(nor, d, k, inside=False)
    wall2 = Wall(nor, d, k, inside=True)
    atom = Atom(np.array([1,0,0]), 1, 'atom')
    pot1 = wall1.calc_potential(atom)
    pot2 = wall2.calc_potential(atom)
    print(pot1, pot2)


def builder11():
    
    atoms = []
    atoms.append(Atom(np.array([0.0,0.0,0.0]), 1, 'a1', np.array([0.0,0.0,1.0])))
    atoms.append(Atom(np.array([0.0,0.0,10.0]), 1, 'a2', np.array([0.0,0.0,-1.0])))
    atoms.append(Atom(np.array([0.0,10.0,0.0]), 1, 'a3', np.array([0.0,-1.0,0.0])))
    a1 = Atom(np.array([10.0,0.0,0.0]), 1, 'a4', np.array([-1.0,0.0,0.0]))
    atoms.append(a1)
    a2 = Atom(np.array([9.0,0.0,0.0]), 1, 'a5', np.array([-1.0,0.0,0.0]))
    atoms.append(a2)
    s=System()
    for a in atoms:
        s.add(a)
    return s, a1, a2


def test4(t=2000):
    
    s, a1, a2 = builder11() 
    time = TimeController()
    integrator = VelocityVerlet()
    euler = Euler()
    p = np.array([VanDerWaals(), Wall(np.array([0.0,0.0,1.0]), 11), Wall(np.array([0.0,0.0,-1.0]), 11),
                  Wall(np.array([0.0,1.0,0.0]), 11), Wall(np.array([0.0,-1.0,0.0]), 11),
                  Wall(np.array([1.0,0.0,0.0]), 11), Wall(np.array([-1.0,0.0,0.0]), 11),
                  Bond(a1, a2, 1.0), Resistance()])

    s = Simulation(time, s, integrator, euler, p)

    writer = PDBWriter('test5.pdb')
    s.compile(t, writer)


def builder5():

    atoms = []
    a1 = Atom(np.array([10.0,0.0,0.0]), 1, 'a4', np.array([-1.0,0.0,0.0]))
    atoms.append(a1)
    a2 = Atom(np.array([9.0,0.0,0.0]), 1, 'a5', np.array([0.0,0.0,0.0]))
    atoms.append(a2)
    s=System()
    for a in atoms:
        s.add(a)
    return s, a1, a2


def test5(t=20000):

    s, a1, a2 = builder5()  
    time = TimeController()
    integrator = VelocityVerlet()
    euler = Euler()
    p = np.array([Bond(a1, a2, 2.3), VanDerWaals()])

    s = Simulation(time, s, integrator, euler, p)

    writer = PDBWriter('test6.pdb')
    s.compile(t, writer)

    
def test(t=20000):

    s = System()
    
    s.add(Atom(np.array([1.0,0.0,0.0]),np.array([1.0,1.0,0.0]),1,'H'))
    s.add(Atom(np.array([1.0,0.0,1.0]),np.array([1.0,0.0,1.0]),1,'H'))
    s.add(Atom(np.array([2.0,0.0,0.0]),np.array([0.0,0.0,1.0]),1,'H'))
    s.add(Atom(np.array([4.0,0.0,0.0]),np.array([0.0,0.0,1.0]),1,'H'))
    s.add(Atom(np.array([3.0,0.0,0.0]),np.array([1.0,1.0,0.0]),16,'O'))
    s.add(Atom(np.array([-5.0,0.0,0.0]),np.array([1.0,0.0,1.0]),2,'He'))
    s.add(Atom(np.array([6.0,1.0,0.0]),np.array([2.0,1.0,1.0]),2,'He'))
    s.add(Atom(np.array([5.0,3.0,2.0]),np.array([0.0,0.0,1.0]),2,'He'))

    pot = [
    VanDerWaals(),
    Wall([1,0,0]),
    Wall([-1,0,0]),
    Wall([0,1,0]),
    Wall([0,-1,0]),
    Wall([0,0,1]),
    Wall([0,0,-1]),
    Bond(s.atoms[0],s.atoms[1]),
    Bond(s.atoms[3],s.atoms[2]),
    ]
    time = TimeController()
    integrator = VelocityVerlet()
    euler = Euler()
    pdbwriter = PDBWriter()
    time_controller = TimeController()

    sim = Simulation(time, s, integrator, euler, pot)
    writer = PDBWriter('test7.pdb')
    sim.compile(t, writer)
    
    '''
    wall = Wall([1,0,0])
    van = VanDerWaals()
    bond = Bond(s.atoms[0],s.atoms[1])
    atom = s.atoms[0]
    print('force: '+str(wall.calc_force(s,atom)))
    print(van.calc_force(s, atom))
    print(bond.calc_force(s, atom))
    '''


def test6(t=10000):

    s = System()
    
    s.add(Atom(np.array([3.0,0.0,0.0]),np.array([1.0,1.0,0.0]),1,'H'))
    s.add(Atom(np.array([3.0,0.0,1.0]),np.array([1.0,0.0,1.0]),1,'H'))
    s.add(Atom(np.array([-2.0,0.0,0.0]),np.array([0.0,0.0,1.0]),1,'H'))
    s.add(Atom(np.array([-3.0,0.0,0.0]),np.array([0.0,0.0,1.0]),1,'H'))
    s.add(Atom(np.array([7.0,0.0,0.0]),np.array([1.0,1.0,0.0]),16,'O'))
    s.add(Atom(np.array([8.0,0.0,0.0]),np.array([1.0,0.0,0.0]),16,'O'))
    s.add(Atom(np.array([-5.0,0.0,0.0]),np.array([1.0,0.0,1.0]),2,'He'))
    s.add(Atom(np.array([6.0,1.0,0.0]),np.array([2.0,1.0,1.0]),2,'He'))
    s.add(Atom(np.array([5.0,3.0,2.0]),np.array([0.0,0.0,1.0]),2,'He'))

    pot = [
    VanDerWaals(),
    Wall([1,0,0]),
    Wall([-1,0,0]),
    Wall([0,1,0]),
    Wall([0,-1,0]),
    Wall([0,0,1]),
    Wall([0,0,-1]),
    Bond(s.atoms[0],s.atoms[1]),
    Bond(s.atoms[3],s.atoms[2]),
    Bond(s.atoms[4],s.atoms[5])
    ]
    time = TimeController()
    integrator = VelocityVerlet()
    euler = Euler()
    pdbwriter = PDBWriter()
    time_controller = TimeController()

    sim = Simulation(time, s, integrator, euler, pot)
    writer = PDBWriter('test9.pdb')
    sim.compile(t, writer)
    

def test(cycles):
    system = System()

    system.add(Atom([1, 0, 0], [1, 1, 0], 1, 'H'))
    system.add(Atom([1, 0, 1], [1, 0, 1], 1, 'H'))
    system.add(Atom([2, 0, 0], [0, 0, 1], 1, 'H'))
    system.add(Atom([4, 0, 0], [0, 0, 1], 1, 'H'))
    system.add(Atom([3, 0, 0], [1, 1, 0], 16, 'O'))
    system.add(Atom([-5, 0, 0], [1, 0, 1], 2, 'He'))
    system.add(Atom([6, 1, 0], [2, 1, 1], 2, 'He'))
    system.add(Atom([5, 3, 2], [0, 0, 1], 2, 'He'))

    potentials = [VanDerWaals(system),
                  Wall([1, 0, 0]),
                  Wall([-1, 0, 0]),
                  Wall([0, 1, 0]),
                  Wall([0, -1, 0]),
                  Wall([0, 0, 1]),
                  Wall([0, 0, -1]),
                  Bond(system.atoms[0], system.atoms[1]),
                  Bond(system.atoms[3], system.atoms[2]),
                  ]
    pdbwriter = PDBWriter()
    time_controller = TimeController()

    simulation = Simulation(system, potentials, VelocityVerlet)
    # simulation.start()
    # simulation.run(nsteps=20000)
    # simulation.pdbwriter.close()
    wal = Wall([1, 0, 0])
    van = VanDerWaals(system)
    bond = Bond(system.atoms[0], system.atoms[1])
    atom = system.atoms[0]
    print(wal.calc_force(atom))
    print(van.calc_force(atom))
    print(bond.calc_force(atom))


test(0)


def test():
    system = System()

    system.add(Atom([1, 0, 0], [1, 1, 0], 1, 'H'))
    system.add(Atom([1, 0, 1], [1, 0, 1], 1, 'H'))
    system.add(Atom([2, 0, 0], [0, 0, 1], 1, 'H'))
    system.add(Atom([4, 0, 0], [0, 0, 1], 1, 'H'))
    system.add(Atom([3, 0, 0], [1, 1, 0], 16, 'O'))
    system.add(Atom([-5, 0, 0], [1, 0, 1], 2, 'He'))
    system.add(Atom([6, 1, 0], [2, 1, 1], 2, 'He'))
    system.add(Atom([5, 3, 2], [0, 0, 1], 2, 'He'))

    potentials = [
        VanDerWaals(system),
        Wall([1, 0, 0]),
        Wall([-1, 0, 0]),
        Wall([0, 1, 0]),
        Wall([0, -1, 0]),
        Wall([0, 0, 1]),
        Wall([0, 0, -1]),
        Bond(system.atoms[0], system.atoms[1]),
        Bond(system.atoms[3], system.atoms[2]),
    ]

    simulation = Simulation(system, potentials, VelocityVerlet)
    simulation.start()
    simulation.run(nsteps=20000)
    simulation.pdbwriter.close()
