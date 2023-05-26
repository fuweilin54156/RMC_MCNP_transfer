# -*- coding:utf-8 -*-
# author: Xiaoyu Guo

from RMC.model.input.base import YMLModelObject as BaseModel
from RMC.model.input.Surface import *
import re


class Tally(BaseModel):

    def __init__(self, meshtally=None, surftally=None, celltally=None, tally_bin=None, unparsed=''):
        if meshtally is None:
            meshtally = []
        if surftally is None:
            surftally = []
        if celltally is None:
            celltally = []
        if tally_bin is None:
            tally_bin = []

        self._meshtally = meshtally
        self._surftally = surftally
        self._celltally = celltally
        self._tally_bin = tally_bin
        # todo 添加其他计数器选项的解析
        self._unparsed = unparsed

    @property
    def meshtally(self):
        return self._meshtally

    @property
    def surftally(self):
        return self._surftally

    @property
    def celltally(self):
        return self._celltally

    @property
    def bin(self):
        return self._tally_bin

    def check(self):
        pass

    # 对tally部分对面有引用部分进行修改
    def pro_tally(self, macrobodies=None):
        segment_pattern = re.compile(r'b(\d+)')  # 类似于segment=b1这种形式
        segment_number = []  # 存储已经解析的bin编号，因为surftally与celltally可能会用同一个segment的bin
        if self.surftally:
            for surftally in self.surftally:
                surftally.surf = Surface.surface_list(surftally.surf, macrobodies)
                if surftally.segment:
                    seg_number = int(segment_pattern.search(surftally.segment).group(1))
                    segment_number.append(seg_number)
                    for bin_model in self.bin:
                        if bin_model.number == seg_number:
                            bin_model.value = Surface.surface_list(bin_model.value, macrobodies)
        if self.celltally:
            for celltally in self.celltally:
                if celltally.segment:
                    seg_number = int(segment_pattern.search(celltally.segment).group(1))
                    if seg_number in segment_number:
                        continue
                    else:
                        for bin_model in self.bin:
                            if bin_model.number == seg_number:
                                bin_model.value = Surface.surface_list(bin_model.value, macrobodies)

    def __str__(self):
        card = 'Tally\n'
        for meshtally in self._meshtally:
            card += str(meshtally)
        for surftally in self._surftally:
            card += str(surftally)
        for celltally in self._celltally:
            card += str(celltally)
        for tally_bin in self._tally_bin:
            card += str(tally_bin)
        card += self._unparsed
        card += '\n\n'
        return card


class MeshTally(BaseModel):
    card_option_types = {
        'MESHTALLY': [int],
        'PARTICLE': [int],
        'TYPE': [int],
        'ENERGY': ['list', float, -1],
        'NORMALIZE': [int],
        'HDF5MESH': [int],
        'GEOMETRY': [int],
        'AXIS': ['list', float, 3],
        'VECTOR': ['list', float, 3],
        'ORIGIN': ['list', float, 3],
        'SCOPE': ['list', int, 3],
        'BOUND': ['list', float, 3],
        'SCOPEX': ['list', int, -1],
        'SCOPEY': ['list', int, -1],
        'SCOPEZ': ['list', int, -1],
        'BOUNDX': ['list', float, -1],
        'BOUNDY': ['list', float, -1],
        'BOUNDZ': ['list', float, -1],
    }

    def __init__(self, tally_id=1, tally_type=1, particle=1, energy=None, normalize=0,
                 hdf5mesh=0, geometry=1,
                 axis=None, vector=None, origin=None,
                 scope=None, bound=None,
                 scopex=None, scopey=None, scopez=None,
                 boundx=None, boundy=None, boundz=None):
        self._id = tally_id
        self._type = tally_type
        self._energy = energy
        self._normalize = normalize
        self._hdf5mesh = hdf5mesh
        self._geometry = geometry
        self._axis = axis
        self._vector = vector
        self._origin = origin
        self._scope = scope
        self._bound = bound
        self._scopex = scopex
        self._scopey = scopey
        self._scopez = scopez
        self._boundx = boundx
        self._boundy = boundy
        self._boundz = boundz
        self._particle = particle

    def add_options(self, options):
        if 'MESHTALLY' in options.keys():
            self._id = options['MESHTALLY']
        if 'TYPE' in options.keys():
            self._type = options['TYPE']
        if 'ENERGY' in options.keys():
            self._energy = options['ENERGY']
        if 'NORMALIZE' in options.keys():
            self._normalize = options['NORMALIZE']
        if 'HDF5MESH' in options.keys():
            self._hdf5mesh = options['HDF5MESH']
        if 'GEOMETRY' in options.keys():
            self._geometry = options['GEOMETRY']
        if 'AXIS' in options.keys():
            self._axis = options['AXIS']
        if 'VECTOR' in options.keys():
            self._vector = options['VECTOR']
        if 'ORIGIN' in options.keys():
            self._origin = options['ORIGIN']
        if 'SCOPE' in options.keys():
            self._scope = options['SCOPE']
        if 'BOUND' in options.keys():
            self._bound = options['BOUND']
        if 'SCOPEX' in options.keys():
            self._scopex = options['SCOPEX']
        if 'SCOPEY' in options.keys():
            self._scopey = options['SCOPEY']
        if 'SCOPEZ' in options.keys():
            self._scopez = options['SCOPEZ']
        if 'BOUNDX' in options.keys():
            self._boundx = options['BOUNDX']
        if 'BOUNDY' in options.keys():
            self._boundy = options['BOUNDY']
        if 'BOUNDZ' in options.keys():
            self._boundz = options['BOUNDZ']
        if 'PARTICLE' in options.keys():
            self._particle = options['PARTICLE']

    @property
    def scope(self):
        return self._scope

    @property
    def bound(self):
        return self._bound

    @property
    def scopex(self):
        return self._scopex

    @property
    def scopey(self):
        return self._scopey

    @property
    def scopez(self):
        return self._scopez

    @property
    def boundx(self):
        return self._boundx

    @property
    def boundy(self):
        return self._boundy

    @property
    def boundz(self):
        return self._boundz

    def check(self):
        assert self._id > 0
        assert self._type > 0

        if self._geometry == 2:
            assert self._axis is not None
            assert self._vector is not None
            assert self._origin is not None

        assert ((self._scope is None) == (self._bound is None) and
                (self._scopex is None) == (self._scopey is None) and
                (self._scopey is None) == (self._scopez is None) and
                (self._scopez is None) == (self._boundz is None) and
                (self._boundz is None) == (self._boundx is None) and
                (self._boundx is None) == (self._boundy is None) and
                (self._scope is None) != (self._scopex is None))

    def __str__(self):
        card = 'Meshtally'
        card += ' ' + str(self._id)
        card += ' type = ' + str(self._type)
        card += ' particle = ' + str(self._particle)
        if self._normalize:
            card += ' normalize = ' + str(self._normalize)
        if self._hdf5mesh:
            card += ' hdf5mesh = ' + str(self._hdf5mesh)

        if self._energy is not None:
            card += ' energy =' + self._list_to_str(self._energy)

        if self._geometry != 1:
            card += ' geometry = ' + str(self._geometry)
            card += ' axis = ' + self._list_to_str(self._axis)
            card += ' vector = ' + self._list_to_str(self._vector)
            card += ' origin = ' + self._list_to_str(self._origin)

        if self._scope is not None:
            card += ' scope = ' + self._list_to_str(self._scope)
            card += ' bound = ' + self._list_to_str(self._bound)
        else:
            card += ' scopex = ' + self._list_to_str(self._scopex)
            card += ' scopey = ' + self._list_to_str(self._scopey)
            card += ' scopez = ' + self._list_to_str(self._scopez)
            card += ' boundx = ' + self._list_to_str(self._boundx)
            card += ' boundy = ' + self._list_to_str(self._boundy)
            card += ' boundz = ' + self._list_to_str(self._boundz)

        card += '\n'
        return card

    @staticmethod
    def _list_to_str(var_list):
        string = ''
        for e in var_list:
            string += ' ' + str(e)
        return string


class SurfTally(BaseModel):
    card_option_types = {
        # todo 完善各类选项
        'SURFTALLY': [int],  # id
        'PARTICLE': [int],
        'TYPE': [int],
        'SURF': ['list', float],
        'CELL': ['list', int],
        'FILTER': ['list', int],
        'INTEGRAL': ['list', int],
        'AREA': ['list', float],
        'VECTOR': ['list', float],
        'GAUSSIAN': ['list', float],
        'ENERGY': ['list', float],
        'COSINE': ['list', float],
        'VALUE': [str],
        'SEGMENT': ['list', int],
        'FLAG': [str],
        'SOURCE': [str],
        'NEST': ['list', int],
        'DOSE': ['list', float, -1],
        'MULTIPLIER': ['list', float, -1],
        'ATTENUATOR': ['list', float, -1],
    }

    def __init__(self, tally_id=None, particle=None, tally_type=None, surf=None, cell=None, tally_filter=None,
                 integral=None, area=None, vector=None, gaussian=None, energy=None, cosine=None, value=None,
                 segment=None, flag=None, source=None, nest=None, dose=None, multiplier=None, attenuator=None,
                 unparsed=None):
        self._id = tally_id
        self._particle = particle
        self._type = tally_type
        self._surf = surf
        self._cell = cell
        self._filter = tally_filter
        self._integral = integral
        self._area = area
        self._vector = vector
        self._gaussian = gaussian
        self._energy = energy
        self._cosine = cosine
        self._value = value
        self._segment = segment
        self._flag = flag
        self._source = source
        self._nest = nest
        self._dose = dose
        self._multiplier = multiplier
        self._attenuator = attenuator
        self._unparsed = unparsed

    def check(self):
        assert self._id > 0
        assert self._type > 0

    @property
    def surf(self):
        return self._surf

    @surf.setter
    def surf(self, surf):
        self._surf = surf

    @property
    def segment(self):
        return self._segment

    @segment.setter
    def segment(self, seg):
        self._segment = seg

    def add_options(self, options):
        if 'SURFTALLY' in options.keys():
            self._id = options['SURFTALLY']
        if 'TYPE' in options.keys():
            self._type = options['TYPE']
        if 'ENERGY' in options.keys():
            self._energy = options['ENERGY']
        if 'FILTER' in options.keys():
            self._filter = options['FILTER']
        if 'INTEGRAL' in options.keys():
            self._integral = options['INTEGRAL']
        if 'CELL' in options.keys():
            self._cell = options['CELL']
        if 'SURF' in options.keys():
            self._surf = options['SURF']
        if 'AREA' in options.keys():
            self._area = options['AREA']
        if 'GAUSSIAN' in options.keys():
            self._gaussian = options['GAUSSIAN']
        if 'VALUE' in options.keys():
            self._value = options['VALUE']
        if 'SEGMENT' in options.keys():
            self._segment = options['SEGMENT']
        if 'FLAG' in options.keys():
            self._flag = options['FLAG']
        if 'NEST' in options.keys():
            self._nest = options['NEST']
        if 'VECTOR' in options.keys():
            self._vector = options['VECTOR']
        if 'SOURCE' in options.keys():
            self._source = options['SOURCE']
        if 'COSINE' in options.keys():
            self._cosine = options['COSINE']
        if 'PARTICLE' in options.keys():
            self._particle = options['PARTICLE']
        if 'DOSE' in options.keys():
            self._multiplier = options['DOSE']
        if 'MULTIPLIER' in options.keys():
            self._multiplier = options['MULTIPLIER']
        if 'ATTENUAOTR' in options.keys():
            self._multiplier = options['ATTENUATOR']

    def __str__(self):
        card = 'Surftally'
        card += ' ' + str(self._id)
        card += ' type = ' + str(self._type)
        card += ' surf =' + ' '.join([str(int(x)) for x in self._surf])
        if self._cell is not None:
            card += ' Cell = ' + ' '.join([str(x) for x in self._cell])
        if self._filter is not None:
            card += ' Filter = ' + ' '.join([str(x) for x in self._filter])
        if self._integral is not None:
            card += ' Integral = ' + ' '.join([str(x) for x in self._integral])
        if self._particle is not None:
            card += ' Particle = ' + str(self._particle)
        if self._area is not None:
            card += ' Area = ' + ' '.join([str(x) for x in self._area])
        if self._vector is not None:
            card += ' Vector = ' + ' '.join([str(x) for x in self._vector])
        if self._gaussian is not None:
            card += ' Gaussian = ' + ' '.join([str(x) for x in self._gaussian])
        if self._energy is not None:
            card += ' Energy = ' + ' '.join([str(x) for x in self._energy])
        if self._cosine is not None:
            card += ' Cosine = ' + ' '.join([str(x) for x in self._cosine])
        if self._value is not None:
            card += ' Value = ' + ' '.join([str(x) for x in self._value])
        if self._segment is not None:
            card += ' Segment = ' + ' '.join([str(x) for x in self._segment])
        if self._flag is not None:
            card += ' Flag = ' + ' '.join([str(x) for x in self._flag])
        if self._nest is not None:
            card += ' Nest = ' + ' '.join([str(x) for x in self._nest])
        if self._source is not None:
            card += ' Source = ' + ' '.join([str(x) for x in self._source])
        if self._dose is not None:
            card += ' Dose = ' + ' '.join([str(x) for x in self._dose])
        if self._multiplier is not None:
            card += ' Multiplier = ' + ' '.join([str(x) for x in self._multiplier])
        if self._attenuator is not None:
            card += ' Attenuator = ' + ' '.join([str(x) for x in self._attenuator])

        card += '\n'
        return card


class CellTally(BaseModel):
    card_option_types = {
        # todo 完善各类选项
        'CELLTALLY': [int],  # id
        'TYPE': [int],
        'PARTICLE': [int],
        'ENERGY': ['list', float],
        'TIME': ['list', float],
        'CELL': ['list', int],
        'FILTER': ['list', int],
        'INTEGRAL': ['list', int],
        'VOLUME': ['list', float],
        'GAUSSIAN': ['list', float],
        'VALUE': [str],
        'SEGMENT': ['list', int],
        'FLAG': [str],
        'NEST': ['list', int],
        'REACTIONRATE': [int],
        'SOURCE': [str],
        'DOSE': ['list', float, -1],
        'MULTIPLIER': ['list', float, -1],
        'ATTENUATOR': ['list', float, -1],
    }

    def __init__(self, tally_id=None, tally_type=None, energy=None, tally_filter=None, integral=None, cell=None,
                 time=None, volume=None, gaussian=None, value=None, segment=None, flag=None, nest=None,
                 reactionrate=None, source=None, particle=None, dose=None, multiplier=None, attenuator=None,
                 unparsed=None):
        self._id = tally_id
        self._type = tally_type
        self._energy = energy
        self._filter = tally_filter
        self._integral = integral
        self._cell = cell
        self._time = time
        self._volume = volume
        self._gaussian = gaussian
        self._value = value
        self._segment = segment
        self._flag = flag
        self._nest = nest
        self._reactionrate = reactionrate
        self._source = source
        self._particle = particle
        self._dose = dose
        self._multiplier = multiplier
        self._attenuator = attenuator

        self._unparsed = unparsed

    @property
    def segment(self):
        return self._segment

    def check(self):
        assert self._id > 0
        assert self._type > 0

    def add_options(self, options):
        if 'CELLTALLY' in options.keys():
            self._id = options['CELLTALLY']
        if 'TYPE' in options.keys():
            self._type = options['TYPE']
        if 'PARTICLE' in options.keys():
            self._particle = options['PARTICLE']
        if 'ENERGY' in options.keys():
            self._energy = options['ENERGY']
        if 'FILTER' in options.keys():
            self._filter = options['FILTER']
        if 'INTEGRAL' in options.keys():
            self._integral = options['INTEGRAL']
        if 'CELL' in options.keys():
            self._cell = options['CELL']
        if 'TIME' in options.keys():
            self._time = options['TIME']
        if 'VOLUME' in options.keys():
            self._volume = options['VOLUME']
        if 'GAUSSIAN' in options.keys():
            self._gaussian = options['GAUSSIAN']
        if 'VALUE' in options.keys():
            self._value = options['VALUE']
        if 'SEGMENT' in options.keys():
            self._segment = options['SEGMENT']
        if 'FLAG' in options.keys():
            self._flag = options['FLAG']
        if 'NEST' in options.keys():
            self._nest = options['NEST']
        if 'REACTIONRATE' in options.keys():
            self._reactionrate = options['REACTIONRATE']
        if 'SOURCE' in options.keys():
            self._source = options['SOURCE']
        if 'DOSE' in options.keys():
            self._dose = options['DOSE']
        if 'MULTIPLIER' in options.keys():
            self._multiplier = options['MULTIPLIER']
        if 'ATTENUATOR' in options.keys():
            self._attenuator = options['ATTENUATOR']

    def __str__(self):
        card = 'Celltally'
        card += ' ' + str(self._id)
        card += ' Type = ' + str(self._type)
        card += ' Cell =' + ' '.join([str(x) for x in self._cell])
        if self._energy is not None:
            card += ' Energy = ' + ' '.join([str(x) for x in self._energy])
        if self._time is not None:
            card += ' Time = ' + ' '.join([str(x) for x in self._time])
        if self._filter is not None:
            card += ' Filter = ' + ' '.join([str(x) for x in self._filter])
        if self._integral is not None:
            card += ' Integral = ' + ' '.join([str(x) for x in self._integral])
        if self._volume is not None:
            card += ' Volume = ' + ' '.join([str(x) for x in self._volume])
        if self._gaussian is not None:
            card += ' Gaussian = ' + ' '.join([str(x) for x in self._gaussian])
        if self._value is not None:
            card += ' Card = ' + ' '.join([str(x) for x in self._value])
        if self._segment is not None:
            card += ' Segment = ' + ' '.join([str(x) for x in self._segment])
        if self._particle is not None:
            card += ' Particle = ' + str(self._particle)
        if self._flag is not None:
            card += ' Flag = ' + ' '.join([str(x) for x in self._flag])
        if self._nest is not None:
            card += ' Nest = ' + ' '.join([str(x) for x in self._nest])
        if self._reactionrate is not None:
            card += ' Reactionrate = ' + str(self._reactionrate)
        if self._source is not None:
            card += ' Source = ' + ' '.join([str(x) for x in self._source])
        if self._dose is not None:
            card += ' Dose = ' + ' '.join([str(x) for x in self._dose])
        if self._multiplier is not None:
            card += ' Multiplier = ' + ' '.join([str(x) for x in self._multiplier])
        if self._attenuator is not None:
            card += ' Attenuator = ' + ' '.join([str(x) for x in self._attenuator])
        card += '\n'
        return card


class Bin(BaseModel):
    card_option_types = {
        'BIN': [int],  # id
        'TYPE': [int],
        'VALUE': ['list', float],
        'BOUND': ['list', float]
    }

    def __init__(self, bin_id=None, tally_type=None, value=None, bound=None):
        self._id = bin_id
        self._type = tally_type
        self._value = value
        self._bound = bound

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        self._value = value

    @property
    def number(self):
        return self._id

    def check(self):
        assert self._type == 2
        assert self._id > 0

    def add_options(self, options):
        if 'BIN' in options.keys():
            self._id = options['BIN']
        if 'TYPE' in options.keys():
            self._type = options['TYPE']
        if 'VALUE' in options.keys():
            self._value = options['VALUE']
        if 'BOUND' in options.keys():
            self._bound = options['BOUND']

    def __str__(self):
        card = 'BIN'
        card += ' ' + str(self._id)
        card += ' Type = ' + str(self._type)
        if self._value is not None:
            card += ' Value = '
            if self._type == 1:
                card += ' '.join([str(x) for x in self._value])
            elif self._type == 2:
                card += ' '.join([str(int(x)) for x in self._value])
        if self._bound is not None:
            card += ' Bound = ' + ' '.join([str(x) for x in self._bound])
        card += '\n'
        return card
