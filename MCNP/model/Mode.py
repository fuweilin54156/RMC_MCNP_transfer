# -*- coding:utf-8 -*-
# author: GYH
# date: 2025-06-25

from MCNP.model.base import YMLModelObject as BaseModel


class Mode(BaseModel):
    card_option_types = {
        # todo 完善各类选项
        'MODE':  ['list', str]
    }

    def __init__(self, mode=None, unparsed=''):
        if mode is not None:
            self.mode = mode

        self._unparsed = unparsed

    def __str__(self):
        card = ''
        if self.mode:
            card = 'MODE '+str(self.mode)
        card += '\n'
        if self._unparsed != '':
            card += self._unparsed
        card += '\n\n'
        return card