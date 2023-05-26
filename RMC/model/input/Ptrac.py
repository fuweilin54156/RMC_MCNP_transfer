from RMC.model.input.base import YMLModelObject as BaseModel
from RMC.model.input.Surface import *


class Ptrac(BaseModel):
    yaml_tag = u'!Ptrac'
    card_option_types = {
        'NEU': [int],
        'ELEC': [int],
        'PHO': [int],
        'SRC': [int],
        'BNK': [int],
        'SUR': [int],
        'COL': [int],
        'TER': [int],
        'FILE': [int],
        'WRITE': [int],
        'MAX': [int],
        'BUFFER': [int],
        'MEPH': [int],
        'X': ['list', float],
        'Y': ['list', float],
        'Z': ['list', float],
        'U': ['list', float],
        'V': ['list', float],
        'W': ['list', float],
        'ERG': ['list', float],
        'WGT': ['list', float],
        'VEL': ['list', float],
        'CELL': ['list', int],
        'SURFACE': ['list', float],
        'CELLTALLY': ['list', int],
        'SURFTALLY': ['list', int],
        'MESHTALLY': ['list', int],
        'POINTTALLY': ['list', int],

    }

    def __init__(self, neu=None, elec=None, pho=None, src=None, bnk=None, sur=None, col=None, ter=None, file=None,
                 write=None, max_out=None, buffer=None, meph=None, x=None, y=None, z=None, u=None, v=None, w=None,
                 erg=None, wgt=None, vel=None, cell=None, surface=None, celltally=None, surftally=None, meshtally=None,
                 pointtally=None):
        self._neu = neu
        self._elec = elec
        self._pho = pho
        self._src = src
        self._bnk = bnk
        self._sur = sur
        self._col = col
        self._ter = ter
        self._file = file
        self._write = write
        self._max = max_out
        self._buffer = buffer
        self._meph = meph
        self._x = x
        self._y = y
        self._z = z
        self._u = u
        self._v = v
        self._w = w
        self._erg = erg
        self._wgt = wgt
        self._vel = vel
        self._cell = cell
        self._surface = surface
        self._celltally = celltally
        self._surftally = surftally
        self._meshtally = meshtally
        self._pointtally = pointtally

    def check(self):
        assert self._file > 0

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self, surface):
        self._surface = surface

    def pro_ptrac(self, macrobodies=None):
        if self.surface:
            self.surface = Surface.surface_list(self.surface, macrobodies)

    def add_options(self, options):
        if 'NEU' in options.keys():
            self._neu = options['NEU']
        if 'ELEC' in options.keys():
            self._elec = options['ELEC']
        if 'PHO' in options.keys():
            self._pho = options['PHO']
        if 'SRC' in options.keys():
            self._src = options['SRC']
        if 'BNK' in options.keys():
            self._bnk = options['BNK']
        if 'SUR' in options.keys():
            self._sur = options['SUR']
        if 'COL' in options.keys():
            self._col = options['COL']
        if 'TER' in options.keys():
            self._ter = options['ter']
        if 'FILE' in options.keys():
            self._file = options['FILE']
        if 'WRITE' in options.keys():
            self._write = options['WRITE']
        if 'MAX' in options.keys():
            self._max = options['MAX']
        if 'BUFFER' in options.keys():
            self._buffer = options['BUFFER']
        if 'MEPH' in options.keys():
            self._meph = options['MEPH']
        if 'X' in options.keys():
            self._x = options['X']
        if 'Y' in options.keys():
            self._y = options['Y']
        if 'Z' in options.keys():
            self._z = options['Z']
        if 'U' in options.keys():
            self._u = options['U']
        if 'V' in options.keys():
            self._v = options['V']
        if 'W' in options.keys():
            self._w = options['W']
        if 'ERG' in options.keys():
            self._erg = options['ERG']
        if 'WGT' in options.keys():
            self._wgt = options['WGT']
        if 'VEL' in options.keys():
            self._vel = options['VEL']
        if 'CELL' in options.keys():
            self._cell = options['CELL']
        if 'SURFACE' in options.keys():
            self._surface = options['SURFACE']
        if 'CELLTALLY' in options.keys():
            self._celltally = options['celltally']
        if 'MESHTALLY' in options.keys():
            self._meshtally = options['MESHTALLY']
        if 'SURFTALLY' in options.keys():
            self._surftally = options['SURFTALLY']
        if 'POINTTALLY' in options.keys():
            self._pointtally = options['POINTTALLY']

    def __str__(self):
        card = 'PTRAC'
        if self._neu is not None:
            card += '  NEU = ' + str(self._neu)
        if self._elec is not None:
            card += ' ELEC = ' + str(self._elec)
        if self._pho is not None:
            card += ' PHO = ' + str(self._pho)
        if self._src is not None:
            card += ' SRC = ' + str(self._src)
        if self._bnk is not None:
            card += ' BNK = ' + str(self._bnk)
        if self._sur is not None:
            card += ' SUR = ' + str(self._sur)
        if self._col is not None:
            card += ' COL = ' + str(self._col)
        if self._ter is not None:
            card += ' TER = ' + str(self._ter)
        card += ' FILE = ' + str(self._file)
        if self._write is not None:
            card += ' WRITE = ' + str(self._write)
        if self._max is not None:
            card += ' MAX = ' + str(self._max)
        if self._buffer is not None:
            card += ' BUFFER = ' + str(self._buffer)
        if self._meph is not None:
            card += ' MEPH = ' + str(self._meph)
        if self._x is not None:
            card += ' X = ' + ' '.join([str(x) for x in self._x])
        if self._y is not None:
            card += ' Y = ' + ' '.join([str(x) for x in self._y])
        if self._z is not None:
            card += ' Z = ' + ' '.join([str(x) for x in self._z])
        if self._u is not None:
            card += ' U = ' + ' '.join([str(x) for x in self._u])
        if self._v is not None:
            card += ' V = ' + ' '.join([str(x) for x in self._v])
        if self._w is not None:
            card += ' W = ' + ' '.join([str(x) for x in self._w])
        if self._erg is not None:
            card += ' ERG = ' + ' '.join([str(x) for x in self._erg])
        if self._wgt is not None:
            card += ' WGT = ' + ' '.join([str(x) for x in self._wgt])
        if self._vel is not None:
            card += ' VEL = ' + ' '.join([str(x) for x in self._vel])
        if self._cell is not None:
            card += ' CELL = ' + ' '.join([str(x) for x in self._cell])
        if self._surface is not None:
            card += ' SURFACE = ' + ' '.join([str(int(x)) for x in self._surface])
        if self._celltally is not None:
            card += ' CELLTALLY = ' + ' '.join([str(x) for x in self._celltally])
        if self._surftally is not None:
            card += ' surftally = ' + ' '.join([str(x) for x in self._surftally])
        if self._meshtally is not None:
            card += ' MESHTALLY = ' + ' '.join([str(x) for x in self._meshtally])
        if self._pointtally is not None:
            card += ' POINTTALYY = ' + ' '.join([str(x) for x in self._pointtally])
        card += '\n\n'
        return card
