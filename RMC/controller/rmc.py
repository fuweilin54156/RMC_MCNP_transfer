# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2019-11-27
import re

from RMC.controller.controller import Controller
from RMC.controller.second_source import SecondSource
from RMC.parser.PlainParser import PlainParser
from RMC.controller.refuel import Refuel
from RMC.parser.YMLParser import YMLParser

from RMC.model.input.Include import IncludeMaterial
from RMC.model.util.material import get_negative_mat_file_map

import os
import shutil
import glob
import re
import warnings


class RMCController(Controller):
    COUPLE_INDEX = 0
    BURNUP_INDEX = 1
    CYCLE_INDEX = 2

    class Options:
        def __init__(self, burnup=False, refuel=False,
                     couple=False, pcqs=False, second_source=False):
            self.burnup = burnup
            self.refuel = refuel
            self.couple = couple
            self.pcqs = pcqs
            self.second_source = second_source

    def __init__(self, inp, archive=None, conti=False, conti_inp=""):
        self.inp = inp
        self.model = PlainParser(inp).parsed

        # define calculation options
        burnup_mode = self.model['burnup'] is not None
        refuel_mode = self.model['refuelling'] is not None
        couple_mode = self.model['criticality'] is not None and self.model['criticality'].couple_option
        pcqs_mode = False
        second_source_mode = False
        if burnup_mode:
            if self.model['burnup'].second_source_step:
                second_source_mode = True
        if refuel_mode:
            self.refuel = Refuel(os.path.join(os.path.dirname(self.inp),
                                              self.model['refuelling'].file))
        else:
            self.refuel = None
        self.prepared_to_refuel = False

        # 计算过程会先临界迭代couple_iteration_num次，而后进行燃耗/临界，所以非耦合情况下需要置零
        couple_iteration_num = self.model['criticality'].max_iteration if couple_mode else 0
        burnup_step_num = self.model['burnup'].step_number if burnup_mode else 0
        # 默认有一个循环
        cycle_num = len(self.refuel.plan) + 1 if refuel_mode else 1

        # 所有index都是0-base的，其中：
        # couple位和cycle位都不能达到max_index，
        # burnup位可以达到max_index，对应的是全部燃耗结束之后的临界计算
        self.names = ["couple", "burnup", "cycle"]
        self._indexes = [0, 0, 0]
        self._max_indexes = [couple_iteration_num, burnup_step_num, cycle_num]

        self.archive = archive

        self.last_archive = self.archive
        self.last_inp = self.inp

        self.options = RMCController.Options(burnup=burnup_mode,
                                             refuel=refuel_mode,
                                             couple=couple_mode,
                                             pcqs=pcqs_mode,
                                             second_source=second_source_mode)

        # todo 完善接续计算功能
        self.conti = conti
        # self.inp_prefix为输入卡名字的前缀，一般就是最初始的输入卡的名字（不要加上路径），此时self.inp为接续输入卡的路径
        # 在进行接续处理之后，会将self.inp恢复到初始输入卡的路径，同时更新了燃耗步数、计算历史等信息
        self.conti_inp = conti_inp
        if self.conti and self.conti_inp == "":
            raise ValueError("The inp prefix should be provided when continuing calculating is called.")

    @property
    def current_couple_step(self):
        return self._indexes[RMCController.COUPLE_INDEX]

    @property
    def current_burnup_step(self):
        return self._indexes[RMCController.BURNUP_INDEX]

    @property
    def current_cycle_index(self):
        return self._indexes[RMCController.CYCLE_INDEX]

    @current_couple_step.setter
    def current_couple_step(self, couple_step):
        self._indexes[RMCController.COUPLE_INDEX] = couple_step

    @current_burnup_step.setter
    def current_burnup_step(self, burnup_step):
        self._indexes[RMCController.BURNUP_INDEX] = burnup_step

    @current_cycle_index.setter
    def current_cycle_index(self, cycle_step):
        self._indexes[RMCController.CYCLE_INDEX] = cycle_step

    @property
    def total_couple_number(self):
        return self._max_indexes[RMCController.COUPLE_INDEX]

    @property
    def total_burnup_number(self):
        return self._max_indexes[RMCController.BURNUP_INDEX]

    @property
    def total_cycle_number(self):
        return self._max_indexes[RMCController.CYCLE_INDEX]

    def couple_move_forward(self):
        self._indexes[RMCController.COUPLE_INDEX] += 1

    def burnup_move_forward(self):
        self._indexes[RMCController.BURNUP_INDEX] += 1
        if self.model["refuelling"] is not None:
            self.model["refuelling"].reduce_step()

    def cycle_move_forward(self):
        self._indexes[RMCController.CYCLE_INDEX] += 1

    def update_calc_params(self, current_model):
        if self.model["criticality"].vary_cycle is not None:
            cur_cycle = self.model["criticality"].vary_cycle[self.current_couple_step - 1]
            current_model["criticality"].population[2] = cur_cycle

    def arrive_at_refuelling_step(self):
        if not self.options.refuel:
            return False
        return self.model["refuelling"].size() > 0 and self.model["refuelling"].top_step() == 0

    def prepare_inp(self, new_inp_path):
        cur_inp = self.inp
        if not self.options.burnup:
            return PlainParser(self.inp).parsed, cur_inp
        # 解析输入卡
        if self.current_burnup_step == 0:
            # 解析初始输入卡
            cur_model = PlainParser(self.inp).parsed
            if os.path.exists(new_inp_path + '.State.h5'):
                os.remove(new_inp_path + '.State.h5')
        elif self.current_burnup_step <= self.total_burnup_number:
            # 解析上个燃耗步生成的接续输入卡
            tail = ".FMTinp.step1"
            cur_inp = os.path.join(self.last_archive, os.path.basename(self.last_inp) + tail)
            cur_model = self._copy_inp_related_files(cur_inp, os.path.dirname(self.inp),
                                                     copy_inp=False, new_inp=new_inp_path)
        else:
            # 计算全部结束
            cur_model = None
        return cur_model, cur_inp

    def _copy_inp_related_files(self, inp, to_dir, copy_inp=True, new_inp=None):
        original_dir = os.path.dirname(inp)
        if new_inp is None:
            new_inp = os.path.basename(inp)
        if copy_inp:
            shutil.copy(inp, os.path.join(to_dir, new_inp))

        model = PlainParser(inp).parsed
        # 将负号材料对应的npy文件拷贝过来
        mat_map = get_negative_mat_file_map(model.geometry.get_univ(0))
        mat_file_list = ["mat_{}.npy".format(mat_map[cell_idx]) for cell_idx in mat_map]
        for f in mat_file_list:
            shutil.copy(os.path.join(original_dir, f), to_dir)
        # 把上个燃耗步生成的material文件copy过来
        shutil.copy(os.path.join(original_dir, model['includematerial'].material), to_dir)
        if 0 < self.current_burnup_step < self.total_burnup_number:
            # 把上个燃耗步生成的包含点燃耗核素信息的文件copy过来
            # 当`self.current_burnup_step==self.total_burnup_number`时，为燃耗全部完成后的临界计算过程，不需要copy
            last_state = os.path.join(original_dir, os.path.basename(inp) + '.State.h5')
            if os.path.exists(last_state):
                shutil.copy(last_state, os.path.join(to_dir, new_inp + '.State.h5'))
        return model

    def next_step(self, prompt=None):
        """

        :return: whether next step exists, and the model,
                 "True" means the controller is moved forward one step.
                 "False" means the calculation has arrived the destination and should stop. (returned model is None)
        """
        if prompt is None:
            prompt = {}
        prompt['inp'], prompt['archive_dir'] = self._new_inp_archive()
        cur_model, cur_inp = self.prepare_inp(new_inp_path=prompt["inp"])

        def return_val(terminate, current_model=None, criticality_only=False, simplified_criticality=False):
            if terminate:
                return False, None

            if current_model is None:
                raise ValueError("current_model should be specified when terminate is False.")

            if simplified_criticality:
                current_model['print'] = None
                current_model['criticalitysearch'] = None

            if criticality_only or simplified_criticality:
                current_model['burnup'] = None
                current_model['refuelling'] = None

            if current_model["burnup"] is not None:
                # 检查燃耗栅元
                bottom_univ = cur_model['geometry'].get_univ(uid=0)
                current_model['burnup'].check_burn_cell(universe=bottom_univ,
                                                        init_burn_cell=self.model['burnup'].burn_cell)
                # 在推进燃耗步之前，修改耦合程序读取的功率文件，实现变功率耦合迭代
                cur_power = cur_model['burnup'].power[0]
                power_path = os.path.join(os.path.dirname(self.inp), 'TotalPower.inp')
                with open(power_path, 'w') as f:
                    f.write(str(cur_power))
            return True, current_model

        if cur_model is None:
            return return_val(terminate=True)

        if self.current_couple_step < self.total_couple_number:
            # 情形：有耦合计算，且耦合迭代还没完成，或非耦合情况下的单次计算过程
            # 处理：不论是否有燃耗/换料，都应该推进耦合迭代，且不需要输出中间结果
            self.couple_move_forward()
            self.update_calc_params(current_model=cur_model)
            return return_val(terminate=False, current_model=cur_model, simplified_criticality=True)

        # 记录实际的上次使用的inp和archive，并将self.last_inp和self.last_archive更新到本次计算使用的inp和archive
        last_inp, last_archive = self.last_inp, self.last_archive
        self.last_inp, self.last_archive = prompt['inp'], prompt['archive_dir']

        # 情形：纯临界计算，或临界耦合到达最后一步
        if not self.options.burnup:
            if self.current_couple_step == self.total_couple_number:
                if self.options.refuel:
                    self.current_couple_step = 0
                    # 纯临界换料、耦合+换料
                    if self.arrive_at_refuelling_step():
                        self._handle_refuel(prompt['inp'])
                        self.update_calc_params(current_model=cur_model)
                        return return_val(terminate=False, current_model=cur_model)
                    else:
                        return return_val(terminate=True)
                # 保证下次能够进入到下面的else中，终止计算
                self.current_couple_step = self.total_couple_number + 1
                self.update_calc_params(current_model=cur_model)
                return return_val(terminate=False, current_model=cur_model, simplified_criticality=False)
            else:
                return return_val(terminate=True)

        # 情形：有燃耗，如果有耦合则已经完成迭代
        # 处理：判断是否需要换料，如果不需要，则通过当前燃耗步判断是否需要终止计算
        self.current_couple_step = 0
        if self.arrive_at_refuelling_step():
            # 设定：换料卡中的燃耗步数字应该是下一个循环的第一步燃耗的编号，即，换料发生在这个数字对应的燃耗步之前
            # 处理：1. 在换料之前，对应的燃耗步不进行燃耗，而是计算临界，即当前循环的最后一个燃耗之后的临界计算
            #      2. 需要提前准备好换料后的输入卡，作为下面正常的循环计算的初始文件

            # 2. 换料修改输入卡，用于下一循环初始的计算
            refuel_path = os.path.join(self.archive, 'cycle{}'.format(self.current_cycle_index + 1), 'refuel')
            os.makedirs(refuel_path, exist_ok=True)
            refuel_inp = os.path.join(refuel_path, os.path.basename(cur_inp))
            # 将上一个循环的最后一个燃耗步生成的最后一个FMT文件复制到refuel文件夹，而后进行refuel操作
            refuel_model = self._copy_inp_related_files(cur_inp, refuel_path)
            refuel_model["refuelling"].pop_step()
            self._handle_refuel(refuel_inp, cur_model=refuel_model, cur_inp=refuel_inp)
            # 将last相关的参数设置为换料后的输入卡和archive文件夹，从而在下一个循环中，能够按照正常的计算迭代流程进行
            self.last_inp, self.last_archive = last_inp, refuel_path

            # 1. 换料前最后的一个临界计算
            self.update_calc_params(current_model=cur_model)
            return return_val(terminate=False, current_model=cur_model, criticality_only=True)

        elif self.current_burnup_step < self.total_burnup_number:
            self.burnup_move_forward()
            if last_archive[-6:] == "refuel":
                # 刚刚完成换料，这是新循环的第一个燃耗步，需要把平衡氙相关的一些参数归零
                cur_model['burnup'].del_succession_entry('RATIO')
                cur_model['burnup'].del_succession_entry('POINTBURNUP')
                cur_model['burnup'].del_succession_entry('CELLCMLTVBURNUP')
            self.update_calc_params(current_model=cur_model)
            return return_val(terminate=False, current_model=cur_model)
        elif self.current_burnup_step == self.total_burnup_number:
            # 最后一个燃耗步后面的临界计算
            self.burnup_move_forward()
            self.update_calc_params(current_model=cur_model)
            return return_val(terminate=False, current_model=cur_model, criticality_only=True)
        else:
            # 实际上这里永远无法到达，因为燃耗步在上一个分支中全部计算完成，所以下一次进入本函数时，在前面prepare_inp之后
            #   cur_model为None，直接终止计算
            return return_val(terminate=True)

    def check(self, status_file):
        """
        Check whether the execution finished by inspecting the status file.
        :param status_file: The path to the status file printed by RMC.
        :return:
        """
        status = YMLParser(status_file).parsed['RMC']
        step = 0
        for cycle_entry in status[:-1]:
            cycle = cycle_entry['cycle']
            step += len(cycle) - 1

        last_cycle = status[-1]['cycle']
        if 'output' in last_cycle[-1]:
            step += len(last_cycle) - 1
            inp = last_cycle[-1]['output']
        else:
            raise NotImplementedError('Burnup restart feature currently '
                                      'not supported in Python Package.')

        if self.refuel is not None:
            # todo: correct the save_remove dictionary here.
            # todo: check directory of inp and whether this check is valid.
            output_file = self.refuel.refuel(step, inp, self.model, save_remove={
                "storage": os.path.dirname(status_file),
                "cycle": -1,
            })
            if output_file is not None:
                return [False, output_file]
            else:
                return [True]
        return [True]

    def _new_inp_archive(self, step=None):
        new_inp = self.inp
        new_archive = self.archive
        if step is None:
            current_burnup_step = self.current_burnup_step
        else:
            current_burnup_step = step

        if self.options.refuel:
            new_inp += '.cycle.' + str(self.current_cycle_index + 1)
            new_archive = os.path.join(
                new_archive, 'cycle' + str(self.current_cycle_index + 1))

        if self.options.burnup:
            new_inp += '.burnup.' + str(current_burnup_step)
            new_archive = os.path.join(new_archive, 'burnup' + str(current_burnup_step))

        if self.options.couple:
            new_inp += '.couple.' + str(self.current_couple_step + 1)
            new_archive = os.path.join(new_archive, 'couple' + str(self.current_couple_step + 1))

        if new_inp == self.inp:
            new_inp = self.inp + '.rmc'

        return [new_inp, new_archive]

    def _handle_whole_burnup(self, prompt):
        if self.current_burnup_step <= self.total_burnup_number:
            # 燃耗步未完成阶段
            prompt['inp'] = self.inp + '.burnup'
            prompt['archive_dir'] = os.path.join(self.archive, 'burnup')
            cur_model = PlainParser(self.inp).parsed
            # 处理二次源强计算
            if self.model['burnup'].second_source_step:
                if cur_model['criticality'] is not None and cur_model['fixedsource'] is not None:
                    cur_model['fixedsource'] = None

            # 直接推进完所有的燃耗步
            self.current_burnup_step = self.total_burnup_number + 1
            # 记录当前燃耗步所用的输入卡和存储位置信息，下个燃耗步会用得到
            self.last_archive = prompt['archive_dir']
            self.last_inp = prompt['inp']
            return True, cur_model
        else:
            # 燃耗步全部完成
            if not self.model['burnup'].second_source_step:
                return False, None
            else:
                # 二次中子源计算在全部燃耗完成后进行，每次会在最后删掉一个二次源计算步
                # todo: 目前不支持接续
                step = self.model['burnup'].second_source_step[0]
                # 获得指定燃耗步的的接续输入卡和archive文件夹位置
                inp, archive = self.last_inp, self.last_archive
                # 解析接续输入卡
                second_source_inp = os.path.join(archive, os.path.basename(inp) + '.FMTinp.step' + str(step))
                cur_model = PlainParser(inp=second_source_inp).parsed
                # 根据指定燃耗步生成的Result.h5文件中的二次源信息，生成固定源模式下的通用源
                second_source_result_h5 = os.path.join(archive, os.path.basename(inp) + '.Result.h5')
                second_source = SecondSource(step=step)
                cur_model['externalsource'] = second_source(second_source_result_h5)
                # 更新二次源启动时计算模型
                # 更新固定源计算条件,更新物理模型
                cur_model['fixedsource'] = self.model['fixedsource']
                cur_model['physics'] = self.model['physics']
                # 关闭燃耗计算中的临界、燃耗及输出模块(二次源启动是固定源模式下的中光子混合输运)
                cur_model['criticality'] = None
                cur_model['burnup'] = None
                cur_model['print'] = None
                # 　获得二次源启动的输入卡及输出目录
                prompt['inp'] = self.inp + '.second_source.' + str(step)
                prompt['archive_dir'] = os.path.join(self.archive, 'second_source' + str(step))
                shutil.copy(
                    os.path.join(self.last_archive,
                                 cur_model['includematerial'].material),
                    os.path.dirname(self.inp)
                )
                del self.model['burnup'].second_source_step[0]
                return True, cur_model

    def continuing(self, status_file, prompt):
        if self.conti:
            self.recover_from_conti()
            # 接续计算部分在每次连续计算中只能运行一次
            self.conti = False

        if self.options.couple:
            if self.current_couple_step == 0:
                # 换料后的第一步也是用正弦分布
                if not self.options.burnup or self.current_burnup_step == 0:
                    # 在耦合计算刚刚开始时，为热工程序生成初始的正弦功率分布
                    # RMC要求用于耦合计算的网格计数器在计数器选项卡中的编号必须是1
                    print('Generating sine shaped power distribution for CTF'
                          'as first guess of the results\n', flush=True)
                    # 在第一次计算之前，先清理先前的计算可能残留的文件
                    trash_files = glob.glob(os.path.join(os.path.dirname(prompt['inp']), 'MeshTally*'))
                    for trash in trash_files:
                        os.remove(trash)
                    # 生成功率
                    self.generate_tally_hdf(self.model['tally'].meshtally[0], os.path.dirname(prompt['inp']))
        if self.model["burnup"] is not None and not self.model['burnup'].single_step:
            succeeded, cur_model = self._handle_whole_burnup(prompt)
        else:
            succeeded, cur_model = self.next_step(prompt)
        if succeeded:
            # 生成实际计算所需要的输入卡
            with open(prompt['inp'], 'w') as f:
                f.write(str(cur_model))
            return True
        else:
            return False

        # if status_file is None:
        #     warnings.warn('Status file not specified, continuous running can not be performed.')
        #     return False
        # status = self.check(status_file)
        # if status[0]:
        #     return False
        # else:
        #     # todo: more elegant method to modify the parameters.
        #     property['inp'] = status[1]
        #     property['conti'] = True
        #     return True

    # TODO: debug时处理step问题
    def _handle_refuel(self, output_file, cur_model=None, cur_inp=None):
        """
        如果换料成功，则current_cycle_index会自动加1

        :param output_file: 换料后的模型输出到的文件路径
        :param cur_model: 换料前的模型
        :param cur_inp: 换料前的输入卡，可以获得一些辅助文件（npy/material等）的路径
        :return: 是否换料成功
        """
        if len(self.model['refuelling'].steps) > 0:
            if self.model['refuelling'].steps[0] == 0:
                step, index = self.model['refuelling'].pop_step()
                # 换料
                save_remove = {
                    'cycle': self.current_cycle_index + 1,
                    'storage': os.path.join(self.archive, "storage")
                }
                if cur_model is None:
                    last_inp = os.path.join(
                        self.last_archive,
                        os.path.basename(self.last_inp)
                    )
                    self.refuel.refuel(index=index, inp=last_inp, base_model=self.model,
                                       output=output_file, dump=True, save_remove=save_remove)
                else:
                    inp = cur_inp if cur_inp is not None else self.inp
                    self.refuel.refuel(index=index, inp=inp, model=cur_model, base_model=self.model,
                                       output=output_file, dump=True, save_remove=save_remove)
                self.cycle_move_forward()
                return True
            else:
                raise ValueError("The step can only be 0 when no burnup")
        else:
            # 换料步列表已空，说明全部换料计算结束
            warnings.warn("Refuel: empty refuelling list in the controller!")
            return False

    @staticmethod
    def generate_tally_hdf(meshtally, h5file_dir):
        """基于解析RMC输入卡得到的网格计数器模型，
        生成一个轴向余弦分布的网格计数器结果HDF5文件，提供给第一次CTF计算使用。

        :param meshtally: 解析RMC输入卡的到的网格计数器模型
        :param h5file_dir: HDF5文件的生成目录
        """
        import numpy as np
        from RMC.controller.RMCEnum import TallyType, MeshType

        def fine_bounds(coarse_bounds, bin_number):
            """ 基于粗网（网格边界坐标和细网格数）计算细网

            :param coarse_bounds: 粗网边界，例如[0.0, 1.0, 3.0]表示粗网有三个边界，
                分别是1.0, 2.0, 3.0
            :param bin_number: 细网格数目，例如[4, 8]表示第一个粗网格中有两个细网格，
                第二个粗网格中有8个细网格
            :return: 细网
            """
            _fine_bounds = []
            for coarse_index in range(len(coarse_bounds) - 1):
                bin_size = coarse_bounds[coarse_index + 1] - \
                           coarse_bounds[coarse_index]  # 细网格的尺寸
                for i in range(bin_number[coarse_index]):
                    # 依次添加细网边界坐标
                    _fine_bounds.append(
                        coarse_bounds[coarse_index] +
                        bin_size / bin_number[coarse_index] * i)
            # 补充上最后的边界坐标
            _fine_bounds.append(coarse_bounds[-1])
            return _fine_bounds

        def sine_power(bounds_x, bounds_y, bounds_z):
            """根据各个方向上的边界（细网），计算所有网格的尺寸（体积）。

            :param bounds_x: x方向上的细网边界，由小到大排列
            :param bounds_y: y方向上的细网边界，由小到大排列
            :param bounds_z: z方向上的细网边界，由小到大排列
            """
            x_bin_num = len(bounds_x) - 1  # mesh bin number in x direction
            y_bin_num = len(bounds_y) - 1  # mesh bin number in y direction
            z_bin_num = len(bounds_z) - 1  # mesh bin number in z direction

            z_bin_center = np.ones(z_bin_num)  # mesh bin center in z direction
            for z_idx in range(z_bin_num):
                z_bin_center[z_idx] = \
                    (bounds_z[z_idx] + bounds_z[z_idx + 1]) / 2.0
            z_length = bounds_z[-1] - bounds_z[0]  # length in z direction

            from math import sin, pi
            # volume of all the meshes
            power = np.ones((x_bin_num, y_bin_num, z_bin_num))
            for x_idx in range(x_bin_num):
                for y_idx in range(y_bin_num):
                    for z_idx in range(z_bin_num):
                        # 功率 = 体积 × 功率密度
                        # 功率密度为正弦分布（堆底为0点）
                        power[x_idx, y_idx, z_idx] = \
                            (bounds_x[x_idx + 1] - bounds_x[x_idx]) * \
                            (bounds_y[y_idx + 1] - bounds_y[y_idx]) * \
                            (bounds_z[z_idx + 1] - bounds_z[z_idx]) * \
                            sin((z_bin_center[z_idx] - bounds_z[0]) / z_length
                                * pi)

            return power

        if meshtally.scope is not None:
            # 均匀结构化网格
            x = meshtally.scope[0]
            y = meshtally.scope[1]
            z = meshtally.scope[2]
            fine_bound_x = fine_bounds([meshtally.bound[0],
                                        meshtally.bound[1]], [x])
            fine_bound_y = fine_bounds([meshtally.bound[2],
                                        meshtally.bound[3]], [y])
            fine_bound_z = fine_bounds([meshtally.bound[4],
                                        meshtally.bound[5]], [z])
        else:
            # 非均匀结构化网格
            x = np.sum(np.array(meshtally.scopex))
            y = np.sum(np.array(meshtally.scopey))
            z = np.sum(np.array(meshtally.scopez))
            fine_bound_x = fine_bounds(meshtally.boundx, meshtally.scopex)
            fine_bound_y = fine_bounds(meshtally.boundy, meshtally.scopey)
            fine_bound_z = fine_bounds(meshtally.boundz, meshtally.scopez)

        bound = []
        bound.extend(fine_bound_x)
        bound.extend(fine_bound_y)
        bound.extend(fine_bound_z)

        power = sine_power(fine_bound_x, fine_bound_y, fine_bound_z)

        import h5py

        # 当前版本的耦合计算要求固定使用下列参数：
        # 用于耦合的网格计数器的编号为1
        mesh_tally_id = 1
        # 用于耦合的网格类型为非均匀的结构化网格
        mesh_type = int(MeshType.nonuniform_structured_mesh)
        # 用于耦合的统计量为功率
        tally_type = int(TallyType.type_power)

        # RMC的网格计数器输出的HDF5文件名为MeshTally{id}.h5
        filename = os.path.join(h5file_dir,
                                u'MeshTally{}.h5'.format(mesh_tally_id))
        h5file = h5py.File(filename, 'w')

        geometry = h5file.create_group('Geometry')
        geometry.attrs['MeshType'] = mesh_type
        geometry.create_dataset('BinNumber', data=np.array([x, y, z]))
        geometry.create_dataset('Boundary', data=np.array(bound))

        h5file.create_dataset('Type' + str(tally_type), data=power)
        h5file.close()

        shutil.copyfile(filename, filename + '.previous')

    def recover_from_conti(self):
        cycle_step, last_burnup_step, last_couple_step = parse_last_calc_info(self.inp, self.conti_inp)
        if self.options.burnup:
            self.current_burnup_step = last_burnup_step + 1
        if self.options.couple:
            # TODO: currently TH coupling can not be continued.
            self.current_couple_step = 0
        if self.options.refuel:
            # cycle index存储为0-base的
            # TODO: currently refuelling can not be continued.
            self.current_cycle_index = cycle_step - 1
        self.last_inp = self.conti_inp

        archive_folder = self.archive
        if self.options.refuel:
            archive_folder = os.path.join(archive_folder, 'cycle{}'.format(cycle_step))
        if self.options.burnup:
            archive_folder = os.path.join(archive_folder, 'burnup{}'.format(last_burnup_step))
        if self.options.couple:
            # TODO: currently TH coupling can not be continued, and the number of the last archive coupling folder is
            #       larger than the max iteration number by 1.
            archive_folder = os.path.join(archive_folder,
                                          "couple{}".format(self.model['criticality'].max_iteration + 1))
        self.last_archive = archive_folder
        print("Continuous calculation recovery report:")
        print("Last cycle step: {}".format(cycle_step))
        print("Last burnup step: {}".format(last_burnup_step))
        print("Last couple step: {}".format(last_couple_step))
        print("Last inp: {}".format(self.last_inp))
        print("Last archive: {}".format(self.last_archive))
        print("Current cycle step: {}".format(self.current_cycle_index + 1))
        print("Current burnup step: {}".format(self.current_burnup_step))


class FakeRMC:
    """假RMC，可以基于RMC的输入卡，生成部分假结果。
    用于辅助对python流程控制代码进行快速debug。
    """

    def __init__(self, inp, archive=None):
        """初始化

        :param inp: RMC的输入卡文件名
        """
        self.inp = inp
        self.archive = archive

    def new_name_start_with(self, base):
        """生成以指定字符串开头的文件名，新文件名以数字结尾区分。

        :param base: 文件名的固定开头
        :return: 新文件名
        """
        if not os.path.exists(os.path.join(os.path.dirname(self.inp), base)):
            return base
        new_name = base
        index = 1
        while os.path.exists(os.path.join(os.path.dirname(self.inp),
                                          new_name)):
            new_name = base + '_' + str(index)
            index += 1
        return new_name

    def archive_output(self):
        if os.path.exists(self.archive):
            print('Warning. Existing folder {} removed.'.format(self.archive))
            shutil.rmtree(self.archive)
        os.makedirs(self.archive)

        file_suffix = ['.out', '.Tally', '.material', '.Adjoint', '.burn.power', '.burn.den_tot',
                       '.depth.error', '.State.h5', '.Result.h5']
        for suffix in file_suffix:
            src = self.inp + suffix
            if not os.path.exists(src):
                continue
            dest = os.path.join(self.archive, os.path.basename(self.inp + suffix))
            shutil.move(src, dest)

        base_dir = os.path.dirname(self.inp)

        for file_folder in os.listdir(base_dir):
            if os.path.isfile(os.path.join(base_dir, file_folder)):
                if re.search(r'material|FMTinp.step\d+', file_folder):
                    dest = os.path.join(self.archive, file_folder)
                    shutil.move(os.path.join(base_dir, file_folder), dest)

    # todo: burnup and refuel step forward
    def run(self):
        """通过解析RMC的输入卡，生成各种输出文件。
        输出的文件中绝大多数都是垃圾内容。
        """
        # dependencies
        assert os.path.exists(os.path.join(os.path.dirname(self.inp), 'xsdir'))

        original_model = PlainParser(self.inp).parsed

        # fake normal results
        with open(self.inp + '.out', 'w') as inpout:
            inpout.write('fake rmc out')
        with open(self.inp + '.Tally', 'w') as tally:
            tally.write('fake rmc tally')
        with open(self.inp + '.material', 'w') as material:
            material.write('fake rmc material')
        with open(self.inp + '.Adjoint', 'w') as adjoint:
            adjoint.write('fake rmc adjoint')

        # fake formatted inp for PRINT block
        if original_model['print'] is not None:
            if original_model['print'].inpfile == 1:
                mat_name = self.new_name_start_with('material')
                with open(os.path.join(os.path.dirname(self.inp),
                                       mat_name), 'w') as material:
                    material.write(str(original_model['material']))
                original_model['material'] = None
                original_model['includematerial'] = IncludeMaterial(mat_name)
                tail = ".FMTinp"
                if original_model['burnup'] is not None:
                    tail += ".step0"
                with open(self.inp + tail, 'w') as fmtinp:
                    fmtinp.write(str(original_model))

        # fake burnup results
        if original_model['burnup'] is not None:
            assert os.path.exists(os.path.join(os.path.dirname(self.inp), 'DepthMainLib'))

            with open(self.inp + '.burn.den_tot', 'w') as dentot:
                dentot.write('fake rmc dentot')
            with open(self.inp + '.burn.power', 'w') as burnpower:
                burnpower.write('fake rmc power')
            with open(self.inp + '.depth.error', 'w') as error:
                error.write('fake rmc error')
            with open(self.inp + '.State.h5', 'w') as state:
                state.write('fake rmc state')
            with open(self.inp + '.Result.h5', 'w') as state:
                state.write('fake rmc result')

        if original_model['burnup'] is not None:
            if original_model['print'] is not None:
                step_number = 2 if original_model['burnup'].single_step else original_model['burnup'].step_number + 1

                for i in range(1, step_number):
                    mat_name = self.new_name_start_with('material')
                    with open(os.path.join(os.path.dirname(self.inp), mat_name), 'w') as material:
                        material.write(str(original_model['material']))

                    original_model['material'] = None
                    original_model['includematerial'] = \
                        IncludeMaterial(mat_name)
                    with open(self.inp + f'.FMTinp.step{i}', 'w') as fmtinp:
                        fmtinp.write(str(original_model))

        mat_dict = get_negative_mat_file_map(original_model.geometry.get_univ(0))
        for i in mat_dict:
            f = "mat_{}.npy".format(mat_dict[i])
            import numpy as np
            a = np.load(f)
            np.save(f, a)

        # print("Fake RMC simulation finished.")


def parse_last_calc_info(inp_base, conti_inp):
    """
    extract step info from the names of input files.

    :param inp_base: original input file
    :param conti_inp: continuous input file
    :return:
    """
    inp_prefix = os.path.basename(inp_base)
    conti_inp_base = os.path.basename(conti_inp)

    length = len(inp_prefix)
    if len(conti_inp_base) < length or conti_inp_base[:length] != inp_prefix:
        raise ValueError("parameters conflict: inp file name {} and prefix {} do not match.".format(conti_inp_base,
                                                                                                    inp_prefix))
    info = conti_inp_base[length:]
    pattern = re.compile(r"^(?:\.cycle\.(\d+))?(?:\.burnup\.(\d+))?(?:\.couple\.(\d+))?$")
    matches = pattern.match(info)
    if not matches:
        raise ValueError("illegal parameter inp file name {}".format(info))
    cycle_step = 0
    burnup_step = 0
    couple_step = 1
    if matches[1] is not None and matches[1] != "":
        cycle_step = int(matches[1])
    if matches[2] is not None and matches[2] != "":
        burnup_step = int(matches[2])
    if matches[3] is not None and matches[3] != "":
        couple_step = int(matches[3])
    return cycle_step, burnup_step, couple_step


if __name__ == '__main__':
    controller = RMCController(inp='inp', archive='archive')
