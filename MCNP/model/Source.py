# -*- coding:utf-8 -*-
# author: GYH
# date: 2025-05-26

from MCNP.model.base import YMLModelObject as BaseModel
import re


class ExternalSource(BaseModel):

    def __init__(self, source=None, distributions=None, unparsed=''):
        if source is None:
            source = []
        if distributions is None:
            distributions = []

        self._source = source
        self._distributions = distributions
        # todo 添加其他源描述选项的解析
        self._unparsed = unparsed

    @property
    def source(self):
        return self._source

    @property
    def distributions(self):
        return self._distributions

    @source.setter
    def source(self, sdef):
        self._source = sdef

    @distributions.setter
    def distributions(self, dist):
        self._distributions = dist

    def __str__(self):
        card = ''
        for source in self._source:
            card += str(source)
        for distribution in self._distributions:
            card += str(distribution)
        card += self._unparsed
        card += '\n\n'
        return card


class Source(BaseModel):
    card_option_types = {
            'ERG': ['list', float],
            'CEL': ['list', int],
            'SUR': ['list', int],
            'TME': ['list', int],
            'DIR': ['list', float],
            'VEC': ['list', float],
            'NRM': [int],
            'POS': ['list', float],
            'RAD': ['list', float],
            'EXT': ['list', float],
            'AXS': ['list', float],
            'X': ['list', float],
            'Y': ['list', float],
            'Z': ['list', float],
            'ARA': [float],
            'WGT': [float],
            'EFF': [float],
            'PAR': ['list', int],
            'TR': ['list', int],
        }

    def __init__(self, erg=None, cell=None, sur=None, time=None, dir=None, vec=None,
                 nrm=None, pos=None, rad=None, ext=None, axs=None, x=None, y=None, z=None, ara=None,
                 wgt=None, eff=None, par=None, tr=None):
        self._erg = erg
        self._cell = cell
        self._sur = sur
        self._time = time
        self._dir = dir
        self._vec = vec
        self._nrm = nrm
        self._pos = pos
        self._rad = rad
        self._ext = ext
        self._axs = axs
        self._x = x
        self._y = y
        self._z = z
        self._ara = ara
        self._wgt = wgt
        self._eff = eff
        self._par = par
        self._tr = tr

    def add_options(self, options):
        if 'ERG' in options.keys():
            self._erg = options['ERG']
        if 'CEL' in options.keys():
            self._cell = options['CEL']
        if 'SUR' in options.keys():
            self._sur = options['SUR']
        if 'TME' in options.keys():
            self._time = options['TME']
        if 'DIR' in options.keys():
            self._dir = options['DIR']
        if 'VEC' in options.keys():
            self._dir = options['VEC']
        if 'POS' in options.keys():
            self._pos = options['POS']
        if 'RAD' in options.keys():
            self._rad = options['RAD']
        if 'EXT' in options.keys():
            self._ext = options['EXT']
        if 'AXS' in options.keys():
            self._axs = options['AXS']
        if 'X' in options.keys():
            self._x = options['X']
        if 'Y' in options.keys():
            self._y = options['Y']
        if 'Z' in options.keys():
            self._z = options['Z']
        if 'ARA' in options.keys():
            self._ara = options['ARA']
        if 'WGT' in options.keys():
            self._wgt = options['WGT']
        if 'EFF' in options.keys():
            self._eff = options['EFF']
        if 'PAR' in options.keys():
            self._par = options['PAR']
        if 'TR' in options.keys():
            self._tr = options['TR']

    def check(self):
        pass

    @property
    def pos(self):
        return self._pos

    @property
    def cell(self):
        return self._cell

    def __str__(self):
        card = 'SDEF'
        if self._erg is not None:
            card += ' ERG = ' + ' '.join([str(x) for x in self._erg])
        if self._cell is not None:
            card += ' CEL = ' + ' '.join([str(x) for x in self._cell])
        if self._sur is not None:
            card += ' SUR = ' + ' '.join([str(x) for x in self._sur])
        if self._time is not None:
            card += ' TME = ' + ' '.join([str(x) for x in self._time])
        if self._dir is not None:
            card += ' DIR = ' + ' '.join([str(x) for x in self._dir])
        if self._vec is not None:
            card += ' VEC = ' + ' '.join([str(x) for x in self._vec])
        if self._nrm is not None:
            card += ' NRM = ' + str(self._nrm)
        if self._pos is not None:
            card += ' POS = ' + ' '.join([str(x) for x in self._pos])
        if self._rad is not None:
            card += ' RAD = ' + ' '.join([str(x) for x in self._rad])
        if self._ext is not None:
            card += ' EXT = ' + ' '.join([str(x) for x in self._ext])
        if self._axs is not None:
            card += ' AXS = ' + ' '.join([str(x) for x in self._axs])
        if self._x is not None:
            card += ' X = ' + ' '.join([str(x) for x in self._x])
        if self._y is not None:
            card += ' Y = ' + ' '.join([str(x) for x in self._y])
        if self._z is not None:
            card += ' Z = ' + ' '.join([str(x) for x in self._z])
        if self._ara is not None:
            card += ' ARA = ' + ' '.join([str(x) for x in self._ara])
        if self._wgt is not None:
            card += ' WGT = ' +  str(self._wgt)
        if self._eff is not None:
            card += ' EFF = ' + ' '.join([str(x) for x in self._eff])
        if self._par is not None:
            card += ' PAR = ' + ' '.join([str(x) for x in self._par])
        if self._tr is not None:
            card += ' TR = ' + ' '.join([str(x) for x in self._tr])
        card += '\n'
        return card


class Distribution(BaseModel):
    card_option_types = {
        'SI': ['list', str],
        'SP': ['list', str],
        'SB': ['list', str],
        'DS': ['list', str],
        'SC': ['list', str]
    }

    def __init__(self, SI=None, SP=None, SB=None, DS=None, SC=None):
        self._SI = SI
        self._SP = SP
        self._SB = SB
        self._DS = DS
        self._SC = SC

    @property
    def SI(self):
        return self._SI

    @SI.setter
    def SI(self, si):
        self._SI = si

    @property
    def SP(self):
        return self._SP

    @SP.setter
    def SP(self, sp):
        self._SP = sp

    @property
    def SB(self):
        return self._SB

    @SB.setter
    def SB(self, sb):
        self._SB = sb

    @property
    def DS(self):
        return self._DS

    @DS.setter
    def DS(self, ds):
        self._DS = ds

    @property
    def SC(self):
        return self._SC

    @SC.setter
    def SC(self, sc):
        self._SC = sc

    def add_options(self, options):
        if 'SI' in options.keys():
            self._SI = options['SI']
        if 'SP' in options.keys():
            self._SP = options['SP']
        if 'SB' in options.keys():
            self._SB = options['SB']
        if 'DS' in options.keys():
            self._DS = options['DS']
        if 'SC' in options.keys():
            self._SC = options['SC']

    def __str__(self):
        card = ''
        if self._SI:
            card += 'SI'
            card += ' '.join([str(x) for x in self._SI])
            card += '\n'
        if self._SP:
            card += 'SP'
            card += ' '.join([str(x) for x in self._SP])
            card += '\n'
        if self._SB:
            card += 'SB'
            card += ' '.join([str(x) for x in self._SB])
            card += '\n'
        if self._DS:
            card += 'DS'
            card += ' '.join([str(x) for x in self._DS])
            card += '\n'
        if self._SC:
            card += 'SC'
            card += ' '.join([str(x) for x in self._SC])
            card += '\n'
        return card

