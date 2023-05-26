# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2019-11-16

from RMC.parser.YMLParser import YMLParser
from RMC.parser.PlainParser import PlainParser
from RMC.utilities.assembly import extract_alias
from RMC.utilities.convert import dict2list
from RMC.model.input.refuelling import FixedPattern
from RMC.model.util.material import get_negative_mat_file_map
from RMC.utilities.assembly import AssemblyRecord

import os
import re
import numpy as np
import warnings
import copy
import h5py


class Refuel:
    def __init__(self, refuel_inp):
        self.base_dir = ""
        self.refuel_inp = refuel_inp
        self.plan = {}
        self.poison_universe = {}
        self.guide_tube = {}
        self.poison_rod = {}
        self.model = None
        self._get_refuel_model()

    def _get_refuel_model(self):
        self.model = YMLParser(self.refuel_inp).parsed['refuelling']
        refuel_list = self.model.lists
        for refuel in refuel_list:
            self.plan[refuel.step] = refuel.plan
            if "poison_universe" in refuel.__dict__.keys():
                self.poison_universe[refuel.step] = refuel.poison_universe
                self.guide_tube[refuel.step] = refuel.guide_tube
                self.poison_rod[refuel.step] = refuel.poison_rod

    def save_as(self, file_name):
        self.model.dump(file_name)

    def refuel(self, index, inp, base_model, save_remove, model=None, output=None, dump=False, base_dir=None):
        # this step is not a refuelling step.
        if index not in self.plan:
            warnings.warn("index {} is not included in the refuelling yaml file, please check whether there are some problems".format(index))
            return None

        self.base_dir = base_dir if base_dir is not None else os.path.dirname(inp)

        plan = self.plan[index]
        if model is None:
            if not os.path.isfile(inp):
                raise FileExistsError("File {} does not exists.".format(inp))
            cur_model = PlainParser(inp).parsed
        else:
            cur_model = model

        # 添加拔棒处理
        if index in self.poison_universe.keys():
            if self.poison_universe[index] is not None:
                for universe_id in self.poison_universe[index]:
                    poison_universe = cur_model['geometry'].get_univ(universe_id)
                    poison_univ_fill_array = poison_universe.lattice.fill
                    poison_univ_fill_array[poison_univ_fill_array == self.poison_rod[index]] = self.guide_tube[index]

        # 删除当前模型中的毒物棒universe，避免读取材料出错
        if index in self.poison_rod.keys():
            cur_model['geometry'].delete_univ(self.poison_rod[index])

        for univ_plan in plan:
            univ_id = univ_plan.universe
            univ = cur_model.geometry.get_univ(univ_id)
            mat_list = self._get_related_mat_file_list(univ=univ)

            univ_pattern = univ_plan.pattern
            univ_list = self.model.get_univ_list(univ_pattern)
            univ_fixed = self.model.get_univ_fixed(univ_pattern)
            univ_alias = self.model.get_univ_alias(univ_pattern)

            exchange_strategy = Refuel._extract_strategy(alias=univ_alias, mapping=univ_plan.mapping)
            new_assembly = set([exchange_strategy[key].assem for key in exchange_strategy
                                if exchange_strategy[key].meta == AssemblyRecord.newAssemType])
            pattern_relationship = \
                Refuel._extract_pattern_relationship(geom=cur_model.geometry, univ=univ,
                                                     universes=univ_list, excepts=univ_fixed,
                                                     extra=new_assembly)
            Refuel._do_refuel(univ_plan.position, univ, mat_list, exchange_strategy, cur_model, base_model,
                              pattern=pattern_relationship, excepts=univ_fixed, save_remove=save_remove)

            # save the universe relationship back.
            self.model.set_univ_list(univ_pattern, dict2list(pattern_relationship[2], null=0))

        if output is not None:
            output_file = output
        else:
            output_file = inp + '.refuel_%d.inp' % index
        with open(output_file, 'w') as f:
            f.write(str(cur_model))

        # todo: write out new refuelling.yml file to extend assembly universes,
        #       or add an option in the universe card to define it is an assembly.

        del (self.plan[index])

        if dump:
            self.save_as(self.refuel_inp)

        return output_file

    @staticmethod
    def _extract_strategy(alias, mapping):
        mapping_data = copy.deepcopy(mapping)
        n_row = len(mapping_data)
        n_column = len(mapping_data[0])
        for i in range(n_row):
            for j in range(n_column):
                if mapping_data[i][j] == 0:
                    mapping_data[i][j] = None
                else:
                    mapping_data[i][j] = extract_alias(alias, mapping_data[i][j])

        strategy = {}
        for i in range(n_row):
            for j in range(n_column):
                initial = i * n_column + j
                if mapping_data[i][j] is not None:
                    strategy[initial] = mapping_data[i][j]

        return strategy

    def _get_related_mat_file_list(self, univ):
        # todo: self.base_dir
        mat_list = get_negative_mat_file_map(univ)
        return {cell_id: os.path.join(self.base_dir, 'mat_{}.npy'.format(mat_list[cell_id])) for cell_id in mat_list}

    @staticmethod
    def _do_refuel(pos, univ, mat_list, strategy, cur_model, base_model, save_remove: dict,
                   pattern=None, excepts: FixedPattern = None):
        initial = univ.lattice.fill.copy()
        new_assemblies = {}

        skip_cycle_assemblies = {}

        if save_remove is not None and 'storage' in save_remove:
            skip_cycle_assemblies = Refuel._extract_skip_cycle_assemblies(strategy, save_remove['storage'])

        pattern_class = pattern[0]  # assembly编号 - 与这个assembly除导管外其他相同的assembly列表
        reverse_pattern = pattern[1]  # assembly编号 - 在上一个dict中对应的key值
        assembly_list = pattern[2]  # core的lattice中的universe编号 - 对应的assembly编号/None
        assembly_univ = pattern[3]  # assembly编号 - 对应的core的lattice中的universe编号

        geom = cur_model.geometry
        for idx in strategy:
            original = univ.lattice.fill[idx]
            if strategy[idx].meta == AssemblyRecord.curCycleType:
                substitute = initial[strategy[idx].pos]
                strategy[idx].assem = substitute
            else:
                substitute = strategy[idx].assem
                if strategy[idx].meta == AssemblyRecord.newAssemType:
                    strategy[idx].pos = idx

            try_update = True
            original_pattern = {'default_pin': [], 'fixed_positions': []}
            substitute_pattern = {'default_pin': [], 'fixed_positions': []}

            if excepts is not None:
                original_pattern['default_pin'] = excepts.get_default_pin(idx)
                original_pattern['fixed_positions'] = excepts.get_positions(idx)
                substitute_pattern['default_pin'] = excepts.get_default_pin(strategy[idx].pos)
                substitute_pattern['fixed_positions'] = excepts.get_positions(strategy[idx].pos)
            if excepts is None or (not excepts.is_fixed_assem(idx) and not excepts.is_fixed_assem(strategy[idx].pos)):
                try_update = False

            if try_update and assembly_list[original] is not None:
                original_assem = geom.get_univ(assembly_list[original])
                substitute_assem = geom.get_univ(assembly_list[substitute])

                default_pin_dict = {}
                subs_idx = 0
                for position in substitute_pattern['fixed_positions']:
                    pin = substitute_pattern['default_pin'][subs_idx]
                    if position[0] not in default_pin_dict:
                        default_pin_dict[position[0]] = {}
                    default_pin_dict[position[0]][position[1]] = pin
                    subs_idx += 1

                for position in original_pattern['fixed_positions']:
                    pin = original_assem.lattice.fill[original_assem.lattice.get_idx(position)]
                    if position[0] not in default_pin_dict:
                        default_pin_dict[position[0]] = {}
                    default_pin_dict[position[0]][position[1]] = pin

                default_pin = []
                for position in excepts.overall_positions:
                    if position[0] not in default_pin_dict or position[1] not in default_pin_dict[position[0]]:
                        default_pin.append(substitute_assem.lattice.fill[substitute_assem.lattice.get_idx(position)])
                    else:
                        default_pin.append(default_pin_dict[position[0]][position[1]])

                fuel, tube = original_assem.compare_pattern(substitute_assem, excepts=excepts.overall_positions)
                if not tube:
                    pattern_id = reverse_pattern[substitute_assem.number]
                    assem_exist = False

                    if not default_pin:
                        default_pin = None

                    for assem in pattern_class[pattern_id]:
                        fuel, tube = original_assem.compare_pattern(geom.get_univ(assem),
                                                                    excepts=excepts.overall_positions,
                                                                    default_pin=default_pin)
                        if tube:
                            substitute = assembly_univ[assem]
                            assem_exist = True
                            break
                    if not assem_exist:
                        if default_pin is None:
                            substitute, end_univ = \
                                geom.get_univ(substitute).duplicate(geom, end_univ=substitute_assem,
                                                                    excepts=excepts.overall_positions,
                                                                    template_univ=original_assem)
                        else:
                            substitute, end_univ = \
                                geom.get_univ(substitute).duplicate(geom, end_univ=substitute_assem,
                                                                    excepts=excepts.overall_positions,
                                                                    default_pin=default_pin)
                        substitute = substitute.number
                        end_univ_id = end_univ.number

                        pattern_class[pattern_id].append(end_univ_id)
                        reverse_pattern[end_univ_id] = pattern_id
                        assembly_list[substitute] = end_univ_id
                        assembly_univ[end_univ_id] = substitute

            if strategy[idx].meta == AssemblyRecord.newAssemType:
                new_assemblies[idx] = substitute
            univ.lattice.fill[idx] = substitute

        except_cells = set()
        if excepts is not None:
            univ_idx_list = excepts.get_fixed_assembly_list()
        else:
            univ_idx_list = []
        for univ_idx in univ_idx_list:
            univ_id = univ.lattice.fill[univ_idx]
            assem_id = assembly_list[univ_id]
            lat = geom.get_univ(assem_id).lattice
            # todo: 删掉重复的检查
            for fill_pos in excepts.get_positions(univ_idx):
                fill_idx = lat.get_idx(fill_pos)
                sub_univ = geom.get_univ(lat.fill[fill_idx])
                for cell in sub_univ:
                    if cell.include is not None:
                        except_mats = Refuel._get_related_mat_file_list(univ=cell.include)
                        for key in except_mats:
                            except_cells.add(key)
                    elif cell.material[0] < 0:
                        except_cells.add(cell.number)

        removed_assemblies = {}
        """
        The expanding method in RMC is DFS, thus those cells belonging to the refuelling universe will be
        neighbors, so only the starting and ending index is needed.
        Also, DFS defines the sorting metric.
        """
        for cell_id in mat_list:
            if cell_id in except_cells:
                continue
            mat_file = mat_list[cell_id]
            mat = np.load(mat_file)
            shape = mat.shape
            """
            0. Find the indexes that will be treated, start and end.
            1. Change the index to move assembly.
            2. Delete some indexes of assemblies that are moved out.
            3. Add material entries for new assemblies.
            4. Sort the cells to fulfill the requirement of RMC lattice material input.
            """
            rev_strategy = Refuel._reverse_strategy(strategy)

            for p in pos:
                p = np.array(p)

                start = -1
                end = -1
                pattern = None
                idx = 0
                while idx <= len(mat):
                    if idx == len(mat):
                        if start >= 0:
                            end = idx
                        else:
                            break
                    else:
                        cell_mat = mat[idx, :]
                        [is_match, pattern] = Refuel._check_pattern_match([pattern, p], cell_mat, univ,
                                                                          model=cur_model)
                        if start < 0:
                            if is_match:
                                start = idx
                            idx += 1
                            continue
                        else:
                            if is_match:
                                idx += 1
                                continue
                            else:
                                end = idx
                                pattern = None

                    aim_cells = mat[start:end, :].copy()
                    mat = np.delete(mat, range(start, end), axis=0)
                    univ_depth = len(p) + 1  # the idx in the mat_row of the No. inside the lattice.
                    to_be_deleted = []

                    for aim_idx in range(len(aim_cells)):
                        cell = aim_cells[aim_idx]
                        if cell[univ_depth] - 1 in rev_strategy:
                            # if the assembly is moved, then change the index
                            # note that after that the cell vector may be wrong
                            # but the sequence of the first dimension is correct.
                            cell[univ_depth] = rev_strategy[cell[univ_depth] - 1] + 1
                        else:
                            if cell[univ_depth] not in removed_assemblies:
                                removed_assemblies[cell[univ_depth]] = \
                                    RefuelAssembly(cycle=save_remove['cycle'],
                                                   position=cell[univ_depth] - 1,
                                                   univ=initial[cell[univ_depth] - 1])
                            removed_assemblies[cell[univ_depth]].append_cell_mats(cell_id, cell[0])
                            to_be_deleted.append(aim_idx)

                    modified_cells = list(np.delete(aim_cells, to_be_deleted, axis=0))

                    # Add material entries for new assemblies.
                    # todo: refactor needed.
                    for assem_pos in new_assemblies:
                        univ_id = new_assemblies[assem_pos]
                        new_univ = cur_model.geometry.get_univ(univ_id)
                        count = new_univ.count_cell(cell_id)
                        if count > 0:
                            initial_mat = base_model.geometry.get_cell(cell_id).material[0]
                            mats = np.zeros((count, shape[1]), dtype=int)
                            mats[:, 0] = -initial_mat
                            mats[:, univ_depth] = assem_pos + 1
                            modified_cells.extend(list(mats))

                    # Add material entries for skip-cycle assemblies.
                    for assem_pos in skip_cycle_assemblies:
                        univ_id = skip_cycle_assemblies[assem_pos].univ
                        new_univ = cur_model.geometry.get_univ(univ_id)
                        assem_mat_list = skip_cycle_assemblies[assem_pos].get_cell_mats(cell_id)
                        count = new_univ.count_cell(cell_id)
                        if count == len(assem_mat_list):
                            mats = np.zeros((count, shape[1]), dtype=int)
                            mats[:, 0] = assem_mat_list
                            mats[:, univ_depth] = assem_pos + 1
                            modified_cells.extend(list(mats))

                    sorted_cells = sorted(modified_cells, key=lambda array: (array[univ_depth]))
                    idx = start + len(sorted_cells)
                    if end - start != len(sorted_cells):
                        warnings.warn("Number of cells removed and inserted does not match for cell {}!".format(cell_id))

                    if sorted_cells:
                        mat = np.insert(mat, start, np.array(sorted_cells), axis=0)
                    start = -1
            np.save(mat_file, mat)

        for pos in removed_assemblies:
            removed_assemblies[pos].dump(save_remove['storage'])

    @staticmethod
    def _reverse_strategy(strategy):
        rev_strategy = {}
        for key in strategy:
            if strategy[key].meta == AssemblyRecord.curCycleType:
                if strategy[key].pos in rev_strategy:
                    raise ValueError('Refuelling matrix error, duplicate assemblies.')
                rev_strategy[strategy[key].pos] = key
        return rev_strategy

    @staticmethod
    def _check_pattern_match(pattern, cell_vec, univ, model):
        length = len(pattern[1])
        if pattern[0] is not None:
            is_match = np.all(pattern[0] == cell_vec[1:length + 1])
            return [is_match, pattern[0]]
        else:
            is_match = True
            for i in range(length):
                if pattern[1][i] > 0 and pattern[1][i] != cell_vec[i + 1]:
                    is_match = False
                    break
            if is_match:
                cell = model.geometry.get_cell(cell_vec[length])
                if cell.include == univ:
                    return [True, np.array(cell_vec[1:length + 1])]
                else:
                    return [False, None]
            else:
                return [False, None]

    @staticmethod
    def _extract_pattern_relationship(geom, univ, universes=None, excepts: FixedPattern = None, extra=None):
        """

        :param geom:
        :param univ:
        :param universes:
        :param excepts:
        :param extra:
        :return:
        """
        extra_list = []
        if extra is not None:
            extra_list = list(extra)
        remain_list = {}

        univ_list = sorted(univ.include, key=lambda x: x.number)
        univ_list.extend([geom.get_univ(uid) for uid in extra_list])

        if universes is not None and universes:
            univ_set = set()
            for pair in universes:
                if len(pair) != 2:
                    raise ValueError("The length of the Universe pair should be 2, this pair is {}".format(len(pair)))
                univ_set.add(pair[1])
                if pair[1] == 0:
                    pair[1] = None
                remain_list[pair[0]] = pair[1]

        for univ in univ_list:
            next_uid = univ.find_next_lattice_univ()
            if next_uid is None:
                warnings.warn("Next universe with lattice not found for Universe {}".format(univ.number))
            remain_list[univ.number] = next_uid

        # note: only one assembly can be corresponding to each lattice universe.
        reverse_list = {}
        assembly_list = []
        for key in remain_list:
            assem = remain_list[key]
            if assem is not None:
                reverse_list[assem] = key
                assembly_list.append(assem)
        assembly_class = {}
        reverse_assembly_class = {}
        for assem in assembly_list:
            new_one = True
            for cls in assembly_class:
                except_positions = None
                if excepts is not None:
                    except_positions = excepts.overall_positions
                fuel, tube = geom.get_univ(cls).compare_pattern(geom.get_univ(assem), excepts=except_positions)
                if fuel:
                    assembly_class[cls].append(assem)
                    reverse_assembly_class[assem] = cls
                    new_one = False
                    break
            if new_one:
                assembly_class[assem] = [assem]
                reverse_assembly_class[assem] = assem

        # assembly_class: [assembly: [assemblies]]
        # reverse_assembly_class: [assembly: assembly_class_id]
        # remain_list: [univ: assembly]
        # remain_list: [assembly: univ]
        return [assembly_class, reverse_assembly_class, remain_list, reverse_list]

    @staticmethod
    def _extract_skip_cycle_assemblies(strategy, archive):
        skip_cycle_assem_map = {}
        for idx in strategy:
            if strategy[idx].meta == AssemblyRecord.skipCycleType:
                cycle = -strategy[idx].assem
                assem = RefuelAssembly.read(archive, cycle, strategy[idx].pos)
                strategy[idx].assem = assem.univ
                skip_cycle_assem_map[idx] = assem
        return skip_cycle_assem_map


class RefuelAssembly:
    pattern = re.compile(r"cell_(\d+)")

    def __init__(self, cycle: int, position: int, univ: int):
        self.cycle = cycle
        self.position = position
        self.univ = univ
        self.burnup_cell_mats = {}

    def set_cell_mats(self, cell_id: int, mats: list):
        self.burnup_cell_mats[cell_id] = np.array(mats, dtype=np.int32)

    # todo: currently does not support top-layer split.
    def append_cell_mats(self, cell_id: int, mat: int):
        if cell_id not in self.burnup_cell_mats:
            self.burnup_cell_mats[cell_id] = []
        self.burnup_cell_mats[cell_id].append(mat)

    def get_cell_mats(self, cell_id):
        if cell_id not in self.burnup_cell_mats:
            return np.array([])
        return self.burnup_cell_mats[cell_id]

    def dump(self, output_dir):
        path = RefuelAssembly.build_path(output_dir, self.cycle)
        if not os.path.exists(path):
            os.makedirs(path)
        elif os.path.isfile(path):
            raise EnvironmentError("Assembly Storage Path \"{}\" should not be occupied by a file".format(path))

        output_file = os.path.join(path, "assem_{}.h5".format(self.position))
        f = h5py.File(output_file, 'w')
        f.attrs.create(name="cycle", data=self.cycle, dtype=np.int32)
        f.attrs.create(name="position", data=self.position, dtype=np.int32)
        f.attrs.create(name="universe", data=self.univ, dtype=np.int32)
        for cell_id in self.burnup_cell_mats:
            f.create_dataset("cell_{}".format(cell_id), data=np.array(self.burnup_cell_mats[cell_id]), dtype=np.int32)
        f.close()

    @staticmethod
    def read(archive, cycle, idx):
        input_file = os.path.join(RefuelAssembly.build_path(archive, cycle), "assem_{}.h5".format(idx))
        f = h5py.File(input_file, 'r')
        assem = RefuelAssembly(cycle=cycle, position=idx, univ=f.attrs.get("universe"))
        for ds in f:
            m = RefuelAssembly.pattern.match(ds)
            if m:
                cell_id = int(m.group(1))
            else:
                raise ValueError("Invalid dataset name in file {}: {}".format(input_file, ds))
            assem.set_cell_mats(cell_id, list(f[ds][...]))
        f.close()
        return assem

    @staticmethod
    def build_path(archive, cycle):
        return os.path.join(archive, "cycle{}".format(cycle))
