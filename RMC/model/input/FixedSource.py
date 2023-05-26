# -*- coding:utf-8 -*-
# author: Hao Luo
# date: 2022-3-15

from RMC.model.input.base import YMLModelObject as BaseModel


class FixedSource(BaseModel):
    card_option_types = {
        'PARTICLE': {
            'POPULATION': [int],
            'FISSION': ['list', int, 2]
        },
        'CUTOFF': {
            'MAXLOST': [int],
            'MINWEIGHT': ['list', float, 2],
            'MAXWEIGHT': [float]
        },
        'LOAD_BALANCE': {
            'METHOD': [bool],
            'INTERVAL': [int]
        }
    }

    def __init__(self, particle=None, cutoff=None, load_balance=None):
        self._particle = particle
        self._cutoff = cutoff
        self._load_balance = load_balance

    @property
    def particle(self):
        return self._particle

    @property
    def cutoff(self):
        return self._cutoff

    @property
    def load_balance(self):
        return self._load_balance

    def check(self):
        pass

    def __str__(self):
        card = 'FIXEDSOURCE\n'
        card += 'Particle population = ' + str(self._particle['POPULATION'])
        if 'FISSION' in self._particle:
            card += ' fission = ' + ' '.join([str(x) for x in self._particle['FISSION']])
        card += '\n'
        if self._cutoff is not None:
            card += 'CutOff'
            if 'MAXLOST' in self._cutoff:
                card += ' maxlost = ' + str(self._cutoff['MAXLOST'])
            if 'MINWEIGHT' in self._cutoff:
                card += ' minweight = ' + ' '.join([str(x) for x in self._cutoff['MINWEIGHT']])
            if 'MAXWEIGHT' in self._cutoff:
                card += ' maxweight = ' + str(self._cutoff['MAXWEIGHT'])
            card += '\n'
        if self._load_balance is not None:
            card += 'Load_Balance'
            if 'METHOD' in self._load_balance:
                card += ' method = ' + '{:d}'.format(self._load_balance['METHOD'])
            if 'INTERVAL' in self._load_balance:
                card += ' interval = ' + str(self._load_balance['INTERVAL'])
            card += '\n'
        card += '\n\n'
        return card
