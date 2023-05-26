# -*- coding:utf-8 -*-
# author: Luo Hao
# date: 2022-10-26

import os
from unittest import TestCase

from RMC.controller.platform import get_platform

from pathlib import Path

dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()


def setUpModule():
    os.chdir(dir_path)


def tearDownModule():
    os.chdir(original_path)


class TestPlatform(TestCase):
    preproc = 'cobratf_preproc'
    ctf = 'cobratf'
    rmc2ctf = 'rmc2ctf'
    ctf2rmc = 'ctf2rmc'
    RMC = 'RMC'
    ctf_n_mpi = 177

    n_mpi = 200
    n_threads = 12
    exec = 'workspace/RMC'
    inp = 'workspace/inp'
    out = None
    restart = 'workspace/inp.binary'
    status = 'workspace/status'
    conti = False
    code = 'RMC'

    reference = 'resources/FakePlatforms/reference'

    def commands(self, platform):
        platform = get_platform(name=platform)
        args_preproc = platform.run_command(mpi=1, commands=self.preproc)
        args_ctf = platform.run_command(mpi=self.ctf_n_mpi, commands=self.ctf)
        args_rmc2ctf = platform.run_command(mpi=1, commands=self.rmc2ctf)
        args_ctf2rmc = platform.run_command(mpi=1, commands=self.ctf2rmc)

        args_rmc = platform.run_command(mpi=self.n_mpi, openmp=self.n_threads, commands=self.exec, inp=self.inp,
                                        out=self.out, restart=self.restart, status=self.status, conti=self.conti,
                                        code="RMC")
        return [args_rmc, args_rmc2ctf, args_ctf, args_preproc, args_ctf2rmc]

    def test_platforms(self):
        linux = self.commands(platform='linux')

        tianhe2 = self.commands(platform='tianhe2')
        tianhe3 = self.commands(platform='tianhe3')
        tianhe3_hpc4 = self.commands(platform='tianhe3_hpc4')
        tianhe3_hpc5 = self.commands(platform='tianhe3_hpc5')
        yinhe = self.commands(platform='yinhe')

        bscc = self.commands(platform='bscc')

        commands_platforms = {'linux': linux, 'tianhe2': tianhe2, 'tianhe3': tianhe3, 'tianhe3_hpc4': tianhe3_hpc4,
                              'tianhe3_hpc5': tianhe3_hpc5, 'yinhe': yinhe, 'bscc': bscc}

        commands = ''

        for key, value in commands_platforms.items():
            commands += f"{key}\n"
            for v in value:
                commands += v + '\n'
            commands += '\n'
        with open(self.reference, 'r') as f:
            reference_commands = f.read()
        self.assertEqual(commands, reference_commands)
