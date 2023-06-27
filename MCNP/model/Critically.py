# -*- coding:utf-8 -*-
# author: GYH
# date: 2025-06-25

from MCNP.model.base import YMLModelObject as BaseModel


class Criticality(BaseModel):
    card_option_types = {
        # todo 完善各类选项
        'KCODE':  ['list', float],
        'KSRC': ['list', float],
    }

    def __init__(self, kcode=None, ksrc=None, unparsed=''):
        if kcode is not None:
            self.kcode = kcode
        if ksrc is not None:
            self.ksrc = ksrc

        self._unparsed = unparsed

    def __str__(self):
        card = ''
        if len(self.kcode) > 0:
            card += 'KCODE '
            card += ' '.join([str(x) for x in self.kcode['KCODE']])
        card += '\n'
        if len(self.ksrc) > 0:
            card += 'KSRC '
            card += ' '.join([str(x) for x in self.ksrc['KSRC']])
        card += '\n'
        if self._unparsed != '':
            card += self._unparsed
        card += '\n\n'
        return card