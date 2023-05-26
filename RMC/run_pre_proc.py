# -*- coding:utf-8 -*-
"""
使用示例

>>> import RMC

>>> RMC.run_file_proc(exec='RMC', inp='pwr_pin', n_mpi=2)

>>> j = RMC.JobFileProc(inp='assembly', n_threads=3)
>>> j.run()

>>> commands = "mpirun -n {n_mpi} {exec} {inp} -s {n_threads}"
>>> RMC.run_file_proc(inp='pwr_pin', n_mpi=2, n_threads=2, commands=commands)

"""

import RMC.preproc.MultiPara as MP
import RMC.preproc.ConstNum as CN
import RMC.preproc.VarNum as VN
import RMC.preproc.FomulaCal as FC
import RMC.preproc.RmDefNotes as RN
import os
import shutil
import RMC
import re
import sys


class JobFileProc:
    def __init__(self, **kwargs):
        self._file_name = os.path.abspath(kwargs['inp'])
        self._kwargs = kwargs
        self._archive_time = {}
        if 'run' in kwargs.keys():  # 是否在处理输入卡后执行计算，默认执行
            self._run = kwargs['run']
        else:
            self._run = True
        if 'status' not in kwargs.keys():  # 提供status初值
            self._kwargs['status'] = None
        if 'archive_dir' not in kwargs.keys():  # 提供archive_dir的初值
            self._kwargs['archive_dir'] = os.path.join(os.getcwd(), 'output')

    def run(self):
        multi_process = MP.MultiPara(self._file_name)
        multi_state = multi_process.devide()

        if multi_state.falseformat:
            print('Multi-value variable format error. Exiting...')
        else:
            try:
                if multi_state.processed:
                    for filename in os.listdir(os.path.join(os.getcwd(), 'MultiParainps')):
                        if re.match(r'.*inp\d+', filename) and not re.search(r'\.', filename):
                            print('Processing ' + filename)
                            shutil.copyfile(os.path.join(os.getcwd(), 'MultiParainps', filename),
                                            os.path.join(os.getcwd(), filename))
                            CN.ConstNum(filename).repconstnum()
                            VN.VarNum(filename).repvar()
                            FC.FomulaCal(filename).calfomula()
                            RN.RmDefNotes(filename).rm_notes()
                            os.makedirs(os.path.join(os.getcwd(), 'MultiParainps', filename + '-output'))
                            if self._run:
                                print(' Start calculating...')
                                self._kwargs.pop('run', 1)
                                self._kwargs['inp'] = os.path.abspath(filename)
                                self._kwargs['archive_dir'] = os.path.join(os.getcwd(), 'MultiParainps',
                                                                           filename + '-output')
                                try:
                                    RMC.run(commands=None, **self._kwargs)
                                except Exception as e:
                                    print(str(e) + ' Exe running error. Exiting...')
                                    sys.exit(1)
                            shutil.copyfile(os.path.join(os.getcwd(), filename),
                                            os.path.join(os.getcwd(), 'MultiParainps', filename + '-output', filename))
                            os.remove(filename)
                else:
                    process_filename = self._file_name + '-cal'
                    shutil.copyfile(self._file_name, process_filename)
                    CN.ConstNum(process_filename).repconstnum()
                    VN.VarNum(process_filename).repvar()
                    FC.FomulaCal(process_filename).calfomula()
                    RN.RmDefNotes(process_filename).rm_notes()
                    self._kwargs['inp'] = process_filename
                    self._kwargs['archive_dir'] = os.path.join(os.getcwd(), process_filename + '-output')
                    if self._run:
                        print(' Start calculating...')
                        self._kwargs.pop('run', 1)
                        try:
                            RMC.run(commands=None, **self._kwargs)
                        except Exception as e:
                            print(str(e) + ' Exe running error. Exiting...')
                            sys.exit(1)
            except Exception as e:
                print(str(e) + ' Format error. Exiting...')
                sys.exit(1)


def run_file_proc(**kwargs):
    job1 = JobFileProc(**kwargs)
    job1.run()


def _argparse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--exec", default="RMC")
    parser.add_argument("--inp", default="inp")
    parser.add_argument("--out", default=None)
    parser.add_argument("--restart", default=None)
    parser.add_argument("--n_mpi", default=None, type=int)
    parser.add_argument("--n_threads", default=None, type=int)
    parser.add_argument("--cwd", default=None)
    parser.add_argument("--print_screen", default=True, type=bool)
    parser.add_argument("--status", default=None)
    parser.add_argument("--conti", default=False, type=bool)
    parser.add_argument("--archive_dir", default=None)
    parser.add_argument("--run", default=True, type=bool)
    parser.add_argument("--proc_per_node", default=None)
    args = parser.parse_args()
    commands = {}
    commands['exec'] = args.exec
    commands['inp'] = args.inp
    commands['out'] = args.out
    commands['restart'] = args.restart
    commands['n_mpi'] = args.n_mpi
    commands['n_threads'] = args.n_threads
    commands['cwd'] = args.cwd
    commands['print_screen'] = args.print_screen
    commands['status'] = args.status
    commands['conti'] = args.conti
    commands['archive_dir'] = args.archive_dir
    commands['run'] = args.run
    commands['proc_per_node'] = args.proc_per_node
    return commands


if __name__ == "__main__":
    commands = _argparse()
    Job1 = JobFileProc(**commands)
    Job1.run()
