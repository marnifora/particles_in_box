import numpy as np
import abc


class Potential(abc.ABC):
    @abc.abstractmethod
    def calc_energy(self, atom1):
        pass

    @abc.abstractmethod
    def calc_force(self, atom1):
        pass


class Wall(Potential):

    def __init__(self, norm_vec, k=2, d=10):
        self.norm_vec = np.array(norm_vec)
        self.k = k
        self.d = d

    def calc_energy(self, atom0):
        r = np.dot(atom0.coords, self.norm_vec) - self.d
        if r < 0:
            return ((-r) ** self.k) * self.norm_vec

        return np.zeros(3)

    def calc_force(self, atom0):
        r = np.dot(atom0.coords, self.norm_vec) - self.d
        if r > 0:
            return (-self.k * (r ** (self.k - 1))) * self.norm_vec

        return np.zeros(3)


class VanDerWaals(Potential):
    def __init__(self, system):
        self.system = system

    def calc_energy(self, atom0):
        potential = np.array([0, 0, 0])
        for atom in self.system.atoms:
            if atom is not atom0:
                r, vec = get_distance(atom0, atom)
                potential += (1 / (r ** 12) - 1 / (r ** 6)) * vec
        return potential

    def calc_force(self, atom0):
        force = np.array([0, 0, 0], dtype='float64')
        for atom in self.system.atoms:
            if atom is not atom0:
                r, vec = get_distance(atom0, atom)
                force += -(6 * (r ** 6 - 2) / r ** 13) * vec
        return force


class Bond(Potential):

    def __init__(self, atom1, atom2, k=1, r0=2.3):
        self.r0 = r0
        self.k = k
        self.atom1 = atom1
        self.atom2 = atom2

    def calc_energy(self, atom0):
        if atom0 is self.atom1:
            r, vec = get_distance(atom0, self.atom2)
        elif atom0 is self.atom2:
            r, vec = get_distance(atom0, self.atom1)
        else:
            return np.zeros(3)
        if r == self.r0:
            return np.zeros(3)
        return self.k * (r - self.r0) ** 2 * vec

    def calc_force(self, atom0):
        if atom0 is self.atom1:
            r, vec = get_distance(atom0, self.atom2)
        elif atom0 is self.atom2:
            r, vec = get_distance(atom0, self.atom1)
        else:
            return np.zeros(3)
        if r == self.r0:
            return np.zeros(3)
        return -self.k * 2 * (r - self.r0) * vec


def get_distance(atom1, atom2):
    AB = atom1.coords - atom2.coords
    norm = np.linalg.norm(AB)
    if norm == 0:
        raise ValueError
    return norm, AB/norm


class LinearPotential(Potential):

    def __init__(self, cx=0, cy=0, cz=0):
        self.cs = np.array([cx, cy, cz])
        
    def calc_force(self, atom):
        return 0

    def calc_energy(self, atom):
        return self.cs[0]*atom.coords[0]+self.cs[1]*atom.coords[1]+self.cs[2]*atom.coords[2]


class RandomPotential(Potential):
    
    def __init__(self):
        
        self.force = None
        
    def calc_force(self, atom):
    
        self.force = np.random.choice(np.arange(0, 1, 0.1), 3)
        return self.force

    def calc_energy(self, atom):
        pass


class Resistance(Potential):

    def __init__(self, f=0.008):
    
        self.f = f
        
    def calc_force(self, atom):
    
        v = np.sqrt(np.sum([i**2 for i in atom.velocity]))
        e = atom.velocity/v
        force = -1*e*self.f*atom.velocity**2
        return force

    def calc_energy(self, atom):
        pass
        
    def get_name(self):
        
        return 'Resistance'


#potencjał wiązania przy inicjalizacji przyjmuje dwa atomy i odległość między nimi, przy liczeniu wartości tego
# potencjału: sprawdza czy odległość pomiędzy zadanymi atomami jest równa zadanej, jeśli nie to liczy potencjał
# zależny od różnicy odległości pomiędzy atomami