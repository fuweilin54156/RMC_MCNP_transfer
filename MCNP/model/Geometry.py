# -*- coding:utf-8 -*-
# author: Shen PF
# date: 2021-07-20
import os

from MCNP.model.base import YMLModelObject as BaseModel
import numpy as np


class Transformation(BaseModel):
    yaml_tag = u'!transformation'

    def __init__(self, paras=None, move=None, rotate=None, num=None, angle=None):
        if move is not None:
            self.move = np.array(move)
        else:
            self.move = None
        if rotate is not None:
            self.rotate = np.array(rotate)
        else:
            self.rotate = None
        self.num = num
        self.paras = paras
        self.angle = angle

    def check(self):
        assert self.move.shape == tuple([3]) and self.rotate.shape == tuple([3, 3])

    def __copy__(self):
        move = None
        if self.move is not None:
            move = np.copy(self.move)
        rotate = None
        if self.rotate is not None:
            rotate = np.copy(self.rotate)
        return Transformation(move=move, rotate=rotate)

    def process(self):
        if len(self.paras) == 1:
            self.num = self.paras[0]
        elif len(self.paras) == 3:
            self.move = np.array(self.paras[0:3])
        elif len(self.paras) == 12:
            self.move = np.array(self.paras[0:3])
            if self.angle is not None:
                self.rotate = np.cos(np.array(self.paras[3:12])*np.pi/180)
            else:
                self.rotate = np.array(self.paras[3:12])
        elif len(self.paras) == 13:
            self.move = np.array(self.paras[0:3])
            if self.angle is not None:
                self.rotate = np.cos(np.array(self.paras[3:12])*np.pi/180)
            else:
                self.rotate = np.array(self.paras[3:12])
            if self.paras[12] == -1:
                self.move = -self.move
        else:
            print("Warning: Please input the total 9 angles or 9 elements in rotate matrix in MCNP input.\n"
                  "The incomplete TR or TRCL card is ignored")


class Cell(BaseModel):
    yaml_tag = u'!cell'
    card_option_types = {
        'FILL': [int],
        'U': [int],
        'LAT': [int],
        'TRCL': ['list', float, -1],
        r'\*TRCL': ['list', float, -1],
        'INNER': [bool],
        'IMP:N': [float],
        'IMP:P': [float],
        'IMP:E': [float]
    }

    def __init__(self, number=-1, bounds='', material=None, density=None, fill=None, inner=False, u=0, lat=None,
                 unparsed=None, impn=None, impp=None, impe=None, trcl=None, likeid=None):
        self.number = number
        self.bounds = bounds
        self.fill = fill
        self.inner = inner
        self.material = material
        self.density = density
        self.universe = u
        self.lat = lat
        self.unparsed = unparsed
        self.impn = impn
        self.impp = impp
        self.impe = impe
        self.trcl = trcl
        self.likeid = likeid

    def check(self):
        assert self.number >= 0

    def add_options(self, options):
        if 'FILL' in options.keys():
            self.fill = options['FILL']
        if 'U' in options.keys():
            self.universe = options['U']
        if 'LAT' in options.keys():
            self.lat = options['LAT']
        if 'IMP:N' in options.keys():
            self.impn = options['IMP:N']
        if 'IMP:P' in options.keys():
            self.impp = options['IMP:P']
        if 'IMP:E' in options.keys():
            self.impe = options['IMP:E']
        if 'TRCL' in options.keys():
            self.trcl = Transformation(paras=options['TRCL'])
        if r'\*TRCL' in options.keys():
            self.trcl = Transformation(paras=options[r'\*TRCL'], angle=1)

    def __str__(self):
        s = '%d %d ' % (self.number, self.material)
        if self.density is not None:
            s += '%f ' % self.density
        else:
            s += '         '
        s += self.bounds + ' '
        if self.universe != 0:
            s += 'u=%d ' % self.universe
        if self.lat is not None:
            s += 'lat=%d ' % self.lat
            s += 'fill=' + ' '.join([str(x) for x in self.fill])
        else:
            if self.fill is not None:
                s += 'fill=%d ' % self.fill
        s += '\n'
        if self.unparsed is not None and self.unparsed != '':
            s += ' Warning: No parserd options : ' + self.unparsed + '\n'
        return s


class Geometry(BaseModel):
    yaml_tag = u'!geometry'

    def __init__(self, universes=None, cell_dict=None, cells=None, unparsed=None):
        if universes is None:
            universes = []
        self.universes = universes
        self.univ_dict = {}
        if cells is None:
            cells = []
        self.cells = cells
        if cell_dict:
            cell_dict = {}
        self.cell_dict = cell_dict
        self.unparsed = unparsed

    def check(self):
        pass

    def add_universe(self, univ):
        self.universes.append(univ)

    def get_univ(self, uid):
        return self.univ_dict[uid]

    def get_cell(self, cid):
        return self.cell_dict[cid]

    def __str__(self):
        s = ''
        for cell in self.cells:
            s += str(cell)
        if self.unparsed is not None and self.unparsed != '':
            s += 'Warning: No parsed cells in geometry block:\n' + self.unparsed
        s += '\n\n'
        return s

    def __iter__(self):
        for univ in self.universes:
            yield univ

    def postprocess(self):
        # scheme1: do the lattice parameter calculations in MCNP block
        pass


class Surfaces(BaseModel):
    yaml_tag = u'!surfaces'

    def __init__(self, surfaces=None, unparsed=''):
        self.surfaces = surfaces
        if self.surfaces is None:
            self.surfaces = []
        self.unparsed = unparsed

    def check(self):
        pass

    def find_surf(self, predicate):
        try:
            return next(idx for idx, n in enumerate(self.surfaces) if predicate(n.number))
        except StopIteration:
            return None

    def update_surface(self, old_surf, new_surf):
        surf_pos = self.find_surf(lambda x: x == old_surf.number)
        if surf_pos is not None:
            assert self.surfaces[surf_pos].type == new_surf.type
            self.surfaces[surf_pos].parameters = new_surf.parameters
        else:
            raise ValueError(
                f'The old surface {old_surf.number} is not in Surfaces of input\n')

    def add_surface(self, surface):
        surf_pos = self.find_surf(lambda x: x == surface.number)
        if surf_pos is None:
            self.surfaces.append(surface)
        else:
            raise ValueError(
                f'The new surface {surface.number} is already in Surfaces of input\n')

    def get_surface(self, surf_num):
        surf_pos = self.find_surf(lambda x: x == surf_num)
        if surf_pos is not None:
            return self.surfaces[surf_pos]
        else:
            raise ValueError(
                f'The surface {surf_num} is not in Surfaces of input')

    def __str__(self):
        s = ''
        for surf in self.surfaces:
            s += str(surf)
        if self.unparsed:
            s += 'Warning: No parsed cards in Surface block: \n' + self.unparsed
        s += '\n\n'
        return s


class Surface(BaseModel):
    yaml_tag = u'!surface'

    surf_type_para = {
        'P': ['list', float, 4],
        'PX': ['list', float, 1],
        'PY': ['list', float, 1],
        'PZ': ['list', float, 1],
        'SO': ['list', float, 1],
        'S': ['list', float, 4],
        'SX': ['list', float, 2],
        'SY': ['list', float, 2],
        'SZ': ['list', float, 2],
        'C/X': ['list', float, 3],
        'C/Y': ['list', float, 3],
        'C/Z': ['list', float, 3],
        'CX': ['list', float, 1],
        'CY': ['list', float, 1],
        'CZ': ['list', float, 1],
        'K/X': ['list', float, 5],
        'K/Y': ['list', float, 5],
        'K/Z': ['list', float, 5],
        'KX': ['list', float, 3],
        'KY': ['list', float, 3],
        'KZ': ['list', float, 3],
        'SQ': ['list', float, 10],
        'GQ': ['list', float, 10],
        'TX': ['list', float, 6],
        'TY': ['list', float, 6],
        'TZ': ['list', float, 6],
        # only for boundary condition
        'BC': [int],
        'PAIR': [int]
    }

    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, tr=None, unparsed=None):
        self.number = number
        self.type = stype
        self.parameters = parameters
        self.boundary = boundary
        self.pair = pair
        self.tr = tr
        self.unparsed = unparsed

    def check(self):
        assert self.number >= 0

    def __str__(self):
        card = str(self.number) + ' ' + self.type + ' '
        surf_para = ' '.join([str(x) for x in self.parameters]) if isinstance(self.parameters, list) else ' ' + str(
            self.parameters)
        card += surf_para
        card += '\n'
        if self.unparsed:
            card += self.unparsed + '\n'
        return card
