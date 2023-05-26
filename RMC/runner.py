# -*- coding:utf-8 -*-
"""Classes/functions for RMC simulation

Examples:

>>> import RMC

>>> RMC.run(exec='RMC', inp='pwr_pin', n_mpi=2)

>>> j = RMC.Job(inp='assembly', n_threads=3)
>>> j.run()

>>> commands = "mpirun -n {n_mpi} {exec} {inp} -s {n_threads}"
>>> RMC.run(inp='pwr_pin', n_mpi=2, n_threads=2, commands=commands)

"""

import time
import subprocess
from numbers import Integral

from RMC.controller.rmc import RMCController
from RMC.controller.rmc import FakeRMC
from RMC.controller.ctf import FakeCTF
from RMC.controller.ctf import FakeCTFPreproc
from RMC.controller.platform import get_platform
from RMC.util.CoupleUtils import power_ave
from RMC.args import parse


class Job(object):
    """A RMC job

    Attributes
    ----------
    exec : str
        Path to RMC executable
    inp : str
        Path to input file
    out : str
        Path to out file
    restart : str
        Path to restart file
    n_mpi : int
        Number of MPI processes
    n_threads : int
        Number of OpenMP threads
    cwd : str
        Path to working directory to run in
    status : str
        Path to status file
    conti : bool
        Whether this run is a continuous one after the previous one
    archive_dir : str
        Path to the directory that contains output / modified files
    platform : str
        Platform of the computer / server. 'linux' or 'tianhe'
    ctf_n_mpi : int
        Number of MPI processes for ctf
    proc_per_node : int
        Number of CPUs for each calculation node of HPC
    """

    def __init__(self, exec='RMC', inp=None, out=None, restart=None,
                 n_mpi=None, n_threads=None, cwd=None, print_screen=True,
                 status=None, conti=False, archive_dir=None, platform='linux',
                 ctf_n_mpi=0, proc_per_node=None, conti_inp="",
                 **kwargs):
        # Initialize class attributes
        self.exec = os.path.abspath(exec)
        self.dir = os.path.dirname(os.path.abspath(inp))
        self.inp = os.path.basename(inp)
        self.out = out
        self.restart = restart
        self.n_mpi = n_mpi
        self.n_threads = n_threads
        self.cwd = cwd
        self.print_screen = print_screen
        self.status = status
        self.conti = conti
        # conti_inp should not be used currently, this parameter is only used to avoid panic for unknown parameters,
        #   the processing of continuing is performed in the controller.
        # TODO: optimize the logic to remove this parameter.
        self.conti_inp = conti_inp

        self.archive_dir = archive_dir
        # 'linux', use mpiexec;
        # 'tianhe', use yhbatch & yhrun;
        # ‘yinhe', use yhbatch & yhrun;
        # 'bscc', use sbatch & srun;
        # bscc 北京超级云计算中心 beijing super cloud computing center
        self.platform = get_platform(name=platform)
        self.ctf_n_mpi = ctf_n_mpi  # 0: shutdown; n: mpi of ctf
        self.cobratf = 'cobratf'
        self.cobratf_preproc = 'cobratf_preproc'
        self.rmc2ctf = 'rmc2ctf'
        self.ctf2rmc = 'ctf2rmc'
        if kwargs:
            raise ValueError('Redundant parameter: ' + str(kwargs.keys()))

    def properties(self):
        properties = ['exec', 'inp', 'out', 'restart', 'n_mpi', 'n_threads', 'status', 'conti']
        return {p: getattr(self, p) for p in properties}

    def __repr__(self):
        return "RMC Job\n {}".format(self.properties())

    @property
    def exec(self):
        return self._exec

    @property
    def inp(self):
        return self._inp

    @property
    def out(self):
        return self._out

    @property
    def restart(self):
        return self._restart

    @property
    def n_mpi(self):
        return self._n_mpi

    @property
    def n_threads(self):
        return self._n_threads

    @property
    def cwd(self):
        return self._cwd

    @property
    def print_screen(self):
        return self._print_screen

    @property
    def status(self):
        return self._status

    @property
    def conti(self):
        return self._conti

    @exec.setter
    def exec(self, exec):
        self._exec = exec

    @inp.setter
    def inp(self, inp):
        self._inp = inp

    @out.setter
    def out(self, out):
        self._out = out

    @restart.setter
    def restart(self, restart):
        self._restart = restart

    @n_mpi.setter
    def n_mpi(self, n_mpi):
        if n_mpi is None or (isinstance(n_mpi, Integral) and n_mpi > 0):
            self._n_mpi = n_mpi
        else:
            raise ValueError(f"Invalid number of MPI processes: {n_mpi}")

    @n_threads.setter
    def n_threads(self, n_threads):
        if n_threads is None or \
                (isinstance(n_threads, Integral) and n_threads > 0):
            self._n_threads = n_threads
        else:
            raise ValueError(f"Invalid number of threads: {n_threads}")

    @cwd.setter
    def cwd(self, cwd):
        self._cwd = cwd

    @print_screen.setter
    def print_screen(self, print_screen):
        self._print_screen = print_screen

    @status.setter
    def status(self, status):
        self._status = status

    @conti.setter
    def conti(self, conti):
        self._conti = conti

    def format_commands(self, commands=None):
        """Generate commands to be run from formatted string"""
        # todo: processing for arguments without parameters, like --continue.
        if commands is None:
            commands = self.platform.run_command(mpi=self.n_mpi, openmp=self.n_threads, commands=self.exec, inp=self.inp,
                                                 out=self.out, restart=self.restart, status=self.status, conti=self.conti, code="RMC")
        return commands.format(**self.properties())

    def _archive_info(self):
        """Obtain the archive information"""
        workspace = self.dir
        all_file_folder = os.listdir(workspace)

        # {'absolute file path': file modification time}
        self._archive_time_stamps = {}
        for file_folder in all_file_folder:
            file_folder = os.path.join(workspace, file_folder)

            # only archive files, not folders
            if os.path.isfile(file_folder):
                self._archive_time_stamps[file_folder] = os.path.getmtime(file_folder)

    def _archive_output(self):
        """Archive output files"""
        import shutil

        all_file_folder = os.listdir(self.dir)

        # create new archive directory
        if os.path.exists(self.archive_dir):
            print('Warning. Existing folder {} removed.'.format(self.archive_dir))
            shutil.rmtree(self.archive_dir)
        os.makedirs(self.archive_dir)

        for file_folder in all_file_folder:
            src = os.path.join(self.dir, file_folder)
            dst = os.path.join(self.archive_dir, file_folder)

            # only archive files, not folders
            if os.path.isfile(src):
                # case 1: new files generated
                if src not in self._archive_time_stamps.keys():
                    shutil.move(src, dst)
                # case 2: files modified
                else:
                    if os.path.getmtime(src) != self._archive_time_stamps[src]:
                        shutil.copy(src, dst)

    def execute(self, args):
        localtime = time.asctime(time.localtime(time.time()))
        print('{} Command to be executed:\n'.format(localtime) + args + '\n')

        # Launch a subprocess
        p = subprocess.Popen(args, cwd=self.cwd, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)

        # Capture and re-print output in real-time
        lines = []
        while True:
            # If run is finished, break loop
            line = p.stdout.readline()
            if not line and p.poll() is not None:
                break

            lines.append(line)
            if self.print_screen:
                # If user requested output, print to screen
                print(line, end='', flush=True)

        # Raise an exception if return status is non-zero
        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, args, ''.join(lines))

    def execute_rmc(self, commands=None):
        args = self.format_commands(commands)
        pwd = os.getcwd()
        os.chdir(self.dir)
        self.execute(args)
        os.chdir(pwd)

    def execute_ctf(self):
        pwd = os.getcwd()
        os.chdir(self.dir)

        preproc = os.path.join(os.getcwd(), self.cobratf_preproc)
        ctf = os.path.join(os.getcwd(), self.cobratf)
        rmc2ctf = os.path.join(os.getcwd(), self.rmc2ctf)
        ctf2rmc = os.path.join(os.getcwd(), self.ctf2rmc)

        assert os.path.exists(preproc), 'ctf preprocessor executable does NOT exists!'
        assert os.path.exists(ctf), 'ctf executable does NOT exists!'
        assert os.path.exists(rmc2ctf), 'interface executable rmc2ctf does NOT exists!'
        assert os.path.exists(ctf2rmc), 'interface executable ctf2rmc does NOT exists!'

        args_preproc = self.platform.run_command(mpi=1, commands=preproc)
        args_ctf = self.platform.run_command(mpi=self.ctf_n_mpi, commands=ctf)
        args_rmc2ctf = self.platform.run_command(mpi=1, commands=rmc2ctf)
        args_ctf2rmc = self.platform.run_command(mpi=1, commands=ctf2rmc)

        self.execute(args_rmc2ctf)
        self.execute(args_preproc)
        self.execute(args_ctf)
        # 改名的原因是，旧版本的CTF在并行时输出pdeck.ctf.h5，而新版本的CTF只输出
        # deck.ctf.h5文件；在接口程序中写死了CTF输出文件的名字为后者
        if os.path.exists(os.path.join(os.getcwd(), 'pdeck.ctf.h5')):
            os.rename(os.path.join(os.getcwd(), 'pdeck.ctf.h5'),
                      os.path.join(os.getcwd(), 'deck.ctf.h5'))
        self.execute(args_ctf2rmc)

        os.chdir(pwd)

    def run(self, commands=None):
        """Run a simulation"""
        self._archive_info()

        pwd = os.getcwd()
        os.chdir(self.dir)

        if self.ctf_n_mpi > 0:
            # 当前限定使用RMC的第一个计数器（必须是符合一定格式的网格计数器）输出的HDF5文件，作为耦合文件
            # 未来有多种热工程序后，上面的if条件会进行调整
            power_ave('MeshTally1.h5')
            self.execute_ctf()
        self.execute_rmc(commands)

        os.chdir(pwd)

        self._archive_output()

    def fake(self):
        """Run a fake simulation"""
        self._archive_info()

        pwd = os.getcwd()
        os.chdir(self.dir)

        if self.ctf_n_mpi > 0:
            power_ave('MeshTally1.h5')
            rmc2ctf = os.path.join(os.getcwd(), 'rmc2ctf')
            args_rmc2ctf = "{}".format(rmc2ctf)
            self.execute(args_rmc2ctf)

            preproc = FakeCTFPreproc(self.inp)
            ctf = FakeCTF(self.inp)

            preproc.run()
            ctf.run()

        rmc = FakeRMC(self.inp)
        rmc.run()

        os.chdir(pwd)

        self._archive_output()


def single_run(commands=None, **kwargs):
    """Run a simulation

    Parameters
    ----------
    commands : str
        Specified string commands to be formatted by the attributes of
        :class:`RMC.Job`
    **kwargs
        Keyword arguments passed to :class:`RMC.Job`

    """
    job = Job(**kwargs)
    job.run(commands=commands)


def run(commands=None, **kwargs):
    """Run a series of RMC simulations.

    Parameters
    ----------
    commands : str
        Specified string commands to be formatted by the attributes of
        :class:`RMC.Job`
    **kwargs
        Keyword arguments passed to :class:`RMC.Job`

    """
    kwargs['inp'] = os.path.abspath(kwargs['inp'])

    controller = RMCController(kwargs['inp'], kwargs['archive_dir'],
                               conti=kwargs['conti'], conti_inp=kwargs['conti_inp'])
    if 'status' not in kwargs:
        print('INFO: No status file specified, simulation terminated...')
        return
    while controller.continuing(kwargs['status'], kwargs):
        single_run(commands=commands, **kwargs)


if __name__ == '__main__':
    import os
    import sys

    os.chdir(sys.path[0])

    archive = os.path.join(os.getcwd(), 'archive')

    args = parse()

    run(inp=args.inp, n_mpi=args.mpi, n_threads=args.omp, exec="workspace/RMC",
        status="workspace/status.txt", conti=args.conti, conti_inp=args.conti_inp, archive_dir=archive,
        ctf_n_mpi=args.assem, platform=args.platform)
