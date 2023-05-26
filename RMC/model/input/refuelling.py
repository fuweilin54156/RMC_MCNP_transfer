# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2019-11-15

from RMC.model.input.base import YMLModelObject as BaseModel
from RMC.utilities.assembly import extract_alias
import yaml
import copy


class RefuelBlock(BaseModel):
    card_option_types = {
        'FILE': [str],
        'STEP': ['list', int, -1],
        'INDEX': ['list', int, -1],
    }

    def __init__(self, file_name=None, steps=None, indexes=None):
        self._file = file_name
        if steps is None:
            self.steps = []
            self.indexes = []
        else:
            if steps != sorted(steps):
                raise ValueError("The refuelling steps should be ordered")
            if len(steps) != len(indexes):
                raise ValueError("The length of refuelling steps and indexes should be equal")
            self.steps = steps
            self.indexes = indexes

        self._current_cycle = 1

    def check(self):
        pass

    def __str__(self):
        card = 'REFUELLING\nFILE ' + self._file
        if len(self.steps) > 0:
            card += '\n'
            card += 'REFUEL step='
            for step in self.steps:
                card += ' ' + str(step)
            card += ' index='
            for index in self.indexes:
                card += ' ' + str(index)
        card += '\n\n'
        return card

    @property
    def file(self):
        return self._file

    def pop_step(self):
        if len(self.steps) == 0:
            raise ValueError("The refueling step vector is empty.")
        else:
            step = self.steps[0]
            index = self.indexes[0]
            self.steps = self.steps[1:]
            self.indexes = self.indexes[1:]
            return step, index

    def reduce_step(self):
        self.steps = [step - 1 for step in self.steps]

    def top_step(self):
        if len(self.steps) == 0:
            raise ValueError("The refueling step vector is empty.")
        else:
            return self.steps[0]

    def size(self):
        return len(self.steps)


class Refuelling(BaseModel):
    yaml_tag = u'!refuelling'

    def __init__(self, lists=None, utilities=None):
        self.lists = lists
        self.utilities = utilities

    def __setstate__(self, state):
        self.lists = state["lists"] if "lists" in state else None
        self.utilities = state["utilities"] if "utilities" in state else None

    def check(self):
        utilities = getattr(self, 'utilities', None)
        if utilities is None:
            setattr(self, 'utilities', {'universes': {}})
        elif 'universes' not in self.utilities:
            self.utilities['universes'] = {}
        pass

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "%s(lists=%r, utilities=%r)" % (self.__class__.__name__, self.lists, self.utilities)

    def __eq__(self, other):
        return repr(self) == repr(other)

    def get_univ_list(self, univ):
        univ_list = self.get_univ_property(univ, 'assemblies')
        if univ_list is None:
            return []
        return univ_list

    def set_univ_list(self, univ, univ_list):
        if univ not in self.utilities['universes']:
            setattr(self, 'utilities', {'universes': {univ: {'assemblies': univ_list}}})
        else:
            self.utilities['universes'][univ]['assemblies'] = univ_list

    def get_univ_fixed(self, univ):
        fixed = self.get_univ_property(univ, 'fixed')
        alias = self.get_univ_alias(univ)
        if fixed is None or alias is None:
            return None
        return FixedPattern(fixed=fixed, alias=alias)

    def get_univ_alias(self, univ):
        # todo: check the shape of the row, column and the mapping.
        alias = self.get_univ_property(univ, 'alias')
        if alias is None:
            return {}
        return alias

    def get_univ_property(self, univ, univ_property):
        if univ not in self.utilities['universes']:
            return None
        elif univ_property in self.utilities['universes'][univ]:
            return self.utilities['universes'][univ][univ_property]
        else:
            return None

    def dump(self, file_name):
        with open(file_name, 'w') as f:
            yaml.dump({'refuelling': self}, f, sort_keys=False, default_flow_style=None)


class Refuel(BaseModel):
    yaml_tag = u'!do_refuel'
    # todo: use __slots__ and types to define __setstate__

    def __init__(self, step=0, plan=None, poison_universe=None, guide_tube=None, poison_rod=None):
        self.step = step
        self.plan = plan
        self.poison_universe = poison_universe
        self.guide_tube = guide_tube
        self.poison_rod = poison_rod

    def check(self):
        pass

    def __setstate__(self, state):
        self.step = state["step"] if "step" in state else 0
        self.plan = state["plan"] if "plan" in state else None
        self.poison_universe = state["poison_universe"] if "poison_universe" in state else None
        self.guide_tube = state["guide_tube"] if "guide_tube" in state else None
        self.poison_rod = state["poison_rod"] if "poison_rod" in state else None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "%s(step=%r, plan=%r, poison_universe=%r, guide_tube=%r, poison_rod=%r)" % \
               (self.__class__.__name__, self.step, self.plan, self.poison_universe, self.guide_tube, self.poison_rod)


class RefuelPlan(BaseModel):
    yaml_tag = u'!refuel_univ'
    alia_name_list = ['column', 'row', 'new']

    def __init__(self, universe=0, position=None, pattern=None, mapping=None):
        self.universe = universe
        self.position = position
        self.mapping = mapping
        self.pattern = pattern

    def __setstate__(self, state):
        self.universe = state["universe"] if "universe" in state else 0
        self.position = state["position"] if "position" in state else None
        self.pattern = state["pattern"] if "pattern" in state else None
        self.mapping = state["mapping"] if "mapping" in state else None

    # todo: check the mapping, number, axis, etc.
    def check(self):
        pass

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "%s(universe=%r, position=%r, pattern=%r, mapping=%r)" % (
            self.__class__.__name__, self.universe, self.position, self.pattern, self.mapping
        )


class FixedPattern:
    def __init__(self, fixed, alias):
        self.assemblies = {}
        self.default_pin = []
        self.positions = []
        self.pos_intervals = []
        self.overall_positions = []
        idx = 0
        pos_start = 0
        for key in fixed:
            for assem_mark in fixed[key]['assemblies']:
                assem_id = extract_alias(alias=alias, univ_coordinate=assem_mark).pos
                if assem_id in self.assemblies:
                    self.assemblies[assem_id].append(idx)
                else:
                    self.assemblies[assem_id] = [idx]
            idx += 1
            if 'default_pin' in fixed[key]:
                self.default_pin.append(fixed[key]['default_pin'])
            else:
                self.default_pin.append(None)
            if 'position' not in fixed[key]:
                raise ValueError('Position parameter should be specified in the fixed block.')
            positions = copy.deepcopy(fixed[key]['position'])
            for pos in positions:
                pos[0] -= 1
                pos[1] -= 1
            self.positions.extend(positions)
            self.pos_intervals.append([pos_start, len(self.positions)])
            pos_start = len(self.positions)
        self._overall_positions()

    def get_default_pin(self, assem_id):
        if self.is_fixed_assem(assem_id):
            result = []
            for idx in self.assemblies[assem_id]:
                interval = self.pos_intervals[idx]
                length = interval[1] - interval[0]
                pin = self.default_pin[idx]
                result.extend([pin] * length)
            return result
        else:
            return []

    def get_positions(self, assem_id):
        if self.is_fixed_assem(assem_id):
            result = []
            for idx in self.assemblies[assem_id]:
                interval = self.pos_intervals[idx]
                result.extend(self.positions[interval[0]:interval[1]])
            return result
        else:
            return []

    def _overall_positions(self):
        pos_dict = {}
        for pos in self.positions:
            if pos[0] not in pos_dict:
                pos_dict[pos[0]] = set()
            pos_dict[pos[0]].add(pos[1])

        self.overall_positions = []
        for key in pos_dict:
            for val in pos_dict[key]:
                self.overall_positions.append([key, val])

    def is_fixed_assem(self, assem_id):
        return assem_id in self.assemblies

    def get_fixed_assembly_list(self):
        return self.assemblies
