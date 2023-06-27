# -*- coding: utf-8 -*-
# author: Xiaoyu Guo

from RMC.model.input.base import YMLModelObject as BaseModel


class Criticality(BaseModel):
    card_option_types = {
        # todo 完善各类选项
        'POWERITER': {
            "KEFF0": [str],  # use str rather than float to keep the original precision
            "POPULATION": ['list', int, 3],
            "BATCHNUM": [int],
        },
        'COUPLE': {
            'MAXITERATION': [int],
            'VARY_CYCLE': ['list', int, -1],
        },
        'INITSRC': {
            'POINT': ['list', float],
            'SLAB': ['list', float],
            'SPH': ['list', float],
            'CLY/X': ['list', float],
            'CYL/Y': ['list', float],
            'CYL/Z': ['list', float],
            'EXTERNALSOURCE': [int],
        }
    }

    def __init__(self, couple=None, power_iter=None, initsrc=None, unparsed=''):
        if couple is not None:
            self._max_iteration = couple["MAXITERATION"] if 'MAXITERATION' in couple else None
            self.vary_cycle = couple["VARY_CYCLE"] if 'VARY_CYCLE' in couple else None
        else:
            self._max_iteration = None
            self.vary_cycle = None
        if power_iter is not None:
            self.keff0 = power_iter["KEFF0"] if 'KEFF0' in power_iter else None
            self.population = power_iter["POPULATION"] if 'POPULATION' in power_iter else None
            self.batch_num = power_iter["BATCHNUM"] if 'BATCHNUM' in power_iter else None
        else:
            raise ValueError("PowerIter option in CRITICALITY card is required.")
        if initsrc is not None:
            self.initsrctype = list(initsrc.keys())[0]
            self.params = initsrc[self.initsrctype]

        self._unparsed = unparsed

        self.couple_option = False
        if self._max_iteration is not None:
            self.couple_option = True

    def check(self):
        if self._max_iteration is not None:
            assert self._max_iteration >= 1
        assert self.population is not None
        if self.vary_cycle is not None:
            assert len(self.vary_cycle) == self.max_iteration

    @property
    def max_iteration(self):
        return self._max_iteration

    @max_iteration.setter
    def max_iteration(self, max_iteration):
        self._max_iteration = max_iteration

    def __str__(self):
        card = 'CRITICALITY\n'
        card += 'PowerIter'
        if self.keff0 is not None:
            card += f' keff0 = {self.keff0}'
        if self.population is not None:
            card += f' population = {int(self.population[0])} {int(self.population[1])} {int(self.population[2])}'
        if self.batch_num is not None:
            card += f' Batchnum = {self.batch_num}'
        card += '\n'
        card += 'InitSrc'
        if self.initsrctype is not None:
            card += f' {self.initsrctype} = '
            if isinstance(self.params,list):
                card += ' '.join([str(x) for x in self.params])
            else:
                card += str(self.params)
        card += '\n'
        if self._unparsed != '':
            card += self._unparsed
        # 注意：max_iteration选项不会输出，因为RMC不需要读取这个选项
        card += '\n\n'
        return card
