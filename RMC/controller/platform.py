# -*- coding:utf-8 -*-
# author: Luo Hao
# date: 2021-5-14

from abc import abstractmethod, ABC
import logging


class Platform:
    def __init__(self, name=None, proc_per_node=None):
        self._name = name
        self._proc_per_node = proc_per_node

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def proc_per_node(self):
        return self._proc_per_node

    @proc_per_node.setter
    def proc_per_node(self, p):
        self._proc_per_node = p

    @abstractmethod
    def run_command(self):
        pass


class Linux(Platform, ABC):
    def __init__(self, name=None, proc_per_node=None):
        super().__init__(name=name, proc_per_node=proc_per_node)

    def run_command(self, **kwargs):
        mpi = kwargs['mpi']
        command = kwargs['commands']
        if 'code' in kwargs.keys() and kwargs['code'] == 'RMC':
            inp = kwargs['inp']
            out = kwargs['out']
            restart = kwargs['restart']
            status = kwargs['status']
            conti = kwargs['conti']

            if 'openmp' in kwargs.keys() and kwargs['openmp'] is not None:
                command += f" -s {kwargs['openmp']} "
            if inp is not None:
                command += inp
            if out is not None:
                command += f" -o {out}"
            if restart is not None:
                command += f" -r {restart}"
            if status is not None:
                command += f"  --status {status}"
            if conti:
                command += " --continue"
        return f"mpiexec -n {mpi} " + command


class TianHeSeries(Platform, ABC):
    def __init__(self, name=None, proc_per_node=None):
        super().__init__(name=name, proc_per_node=proc_per_node)

    def run_command(self, **kwargs):
        command = kwargs['commands']
        if 'code' in kwargs.keys() and kwargs['code'] == 'RMC':
            inp = kwargs['inp']
            out = kwargs['out']
            restart = kwargs['restart']
            status = kwargs['status']
            conti = kwargs['conti']
            if 'openmp' in kwargs.keys() and kwargs['openmp'] is not None:
                command += f" -s {kwargs['openmp']} "
                command = f"-c {kwargs['openmp']} " + command
            if inp is not None:
                command += inp
            if out is not None:
                command += f" -o {out}"
            if restart is not None:
                command += f" -r {restart}"
            if status is not None:
                command += f"  --status {status}"
            if conti:
                command += " --continue"

        n_mpi_omp = kwargs['mpi']
        if 'openmp' in kwargs.keys():
            n_mpi_omp *= kwargs['openmp']
        nodes = n_mpi_omp // self.proc_per_node
        if n_mpi_omp % self.proc_per_node > 0:
            nodes += 1
            print('Warning: node utilization is insufficient. '
                  '{} mpi/threads idle.'.format(nodes * self.proc_per_node - n_mpi_omp))
        command = "yhrun -N " + str(nodes) + f" -n {kwargs['mpi']} " + command
        return command


class BSCC(Platform, ABC):
    def __init__(self, name=None, proc_per_node=None):
        super().__init__(name=name, proc_per_node=proc_per_node)

    def run_command(self, **kwargs):
        command = kwargs['commands']
        if 'code' in kwargs.keys() and kwargs['code'] == 'RMC':
            inp = kwargs['inp']
            out = kwargs['out']
            restart = kwargs['restart']
            status = kwargs['status']
            conti = kwargs['conti']

            if 'openmp' in kwargs.keys() and kwargs['openmp'] is not None:
                command += f" -s {kwargs['openmp']} "
                command = f"-c {kwargs['openmp']} " + command
            if inp is not None:
                command += inp
            if out is not None:
                command += f" -o {out}"
            if restart is not None:
                command += f" -r {restart}"
            if status is not None:
                command += f"  --status {status}"
            if conti:
                command += " --continue"
        return f"srun -n {kwargs['mpi']} " + command


platforms = {"linux": Linux(name="linux"),
             "tianhe2": TianHeSeries(name="tianhe2", proc_per_node=24),
             "tianhe3": TianHeSeries(name="tianhe3", proc_per_node=32),
             "yinhe": TianHeSeries(name="yinhe", proc_per_node=20),
             "tianhe3_hpc4": TianHeSeries(name="tianhe3_hpc4", proc_per_node=36),
             "tianhe3_hpc5": TianHeSeries(name="tianhe3_hpc5", proc_per_node=56),
             "bscc": BSCC(name="bscc", proc_per_node=64),
             "shuguang": Linux(name="shuguang", proc_per_node=12)
             }


def get_platform(name):
    if name in platforms:
        return platforms[name]
    else:
        logging.error(f"Platform {name} cannot be recognized, commands to be run will be set as Linux mode.")
        raise NotImplementedError("Unrecognized calculating platform.")
