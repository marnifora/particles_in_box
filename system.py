from heapq import *
import sys


class System:

    def __init__(self):
        self.atoms = []

    def add(self, atom):
        self.atoms.append(atom)

    def remove(self, atom):
        self.atoms.pop(atom)


class PDBWriter:

    def __init__(self, fname=None):

        if fname is None:
            self.out = sys.stdout
        else:
            self.out = open(fname, 'w')

    def drop_current_state(self, system):

        self.out.write('MODEL\n')
        i = 0
        for a in system.atoms:
            i += 1
            line = 'ATOM  {0:>5} {5:<9}{1}{0:>4}    {2:>8.3f}{3:>8.3f}{4:>8.3f}{1:>6}{5:>18}  \n'.format(
                a.num, 'a', *a.coords, a.name)
            self.out.write(line)
        # for x in range(8):
        #   self.out.write('1234567890')
        self.out.write('ENDMDL\n')

    def close(self):
        self.out.close()


class TimeController:
    def __init__(self):
        self.t = 0
        self.count = 0
        self.event_heap = []

    def add_event(self, time, method, args=[], kwargs={}):
        heappush(self.event_heap, (self.t+time, self.count, method, args, kwargs))
        self.count += 1

    def exec_step(self, step):

        if len(self.event_heap) != 0:
            event = heappop(self.event_heap)
            while event[0] <= self.t:
                event[2](*event[3], **event[4])
                try:
                    event = heappop(self.event_heap)
                except IndexError:
                    event = False
                    break
            if event:
                heappush(self.event_heap, event)
        self.t += step
