# -*- coding:utf-8 -*-
import h5py
import re

from RMC.model.input.base import YMLModelObject as BaseModel
from RMC.model.input.Surface import *


class ExternalSource(BaseModel):
    card_option_types = {
        'H5SOURCEINFO': {
            'FILENAME': [str],
            'GROUPNAME': [str]
        }
    }

    def __init__(self, source=None, surfsrcread=None, distributions=None, h5source=None, unparsed=''):
        if source is None:
            source = []
        if surfsrcread is None:
            surfsrcread = []
        if distributions is None:
            distributions = []

        self._source = source
        self._surfsrcread = surfsrcread
        self._distributions = distributions
        self._h5source = h5source
        # todo 添加其他源描述选项的解析
        self._unparsed = unparsed

    @property
    def source(self):
        return self._source

    @property
    def surfsrcread(self):
        return self._surfsrcread

    @property
    def distributions(self):
        return self._distributions

    @property
    def h5source(self):
        return self._h5source

    @distributions.setter
    def distributions(self, dist):
        self._distributions = dist

    def check(self):
        source_ids = []
        for source in self._source:
            source_ids.append(source._id)
        assert len(set(source_ids)) == len(source_ids), "repeated source id in ExternalSource"

        distribution_ids = []
        for distribution in self._distributions:
            distribution_ids.append(distribution.id)
        assert len(set(distribution_ids)) == len(distribution_ids), "repeated distribution id in ExternalSource"

    def postprocess(self):
        """ 读取HDF5格式的通用源信息"""
        if self._h5source is None:
            return
        source = h5py.File(self._h5source['FILENAME'], 'r')
        group = source['Burnup'][self._h5source['GROUPNAME']]['GammaSource']

        ids = []
        for source in self._source:
            ids.append(source._id)
        for dist in self._distributions:
            ids.append(dist.id)

        id = max(ids)

        for cell in group.keys():
            id += 1
            energy = list(group[cell]['Energy'][...])
            intensity = list(group[cell]['Intensity'][...])
            self._source.append(Source(source_id=id, fraction=1, particle=[2], energy=[f'd{id}'], cell=cell.split()))
            self._distributions.append(Distribution(id=id, type=2, value=energy, probability=intensity))

    def pro_source(self, macrobodies=None):
        if self.source:
            for source in self.source:
                if source.surface:
                    try:
                        source.surface = Surface.surface_list(source.surface, macrobodies)
                    except AttributeError:
                        continue
        if self.surfsrcread:
            if self.surfsrcread[0].newsurf:
                self.surfsrcread[0].newsurf = Surface.surface_list(self.surfsrcread[0].newsurf, macrobodies)
            if self.surfsrcread[0].oldsurf:
                self.surfsrcread[0].oldsurf = Surface.surface_list(self.surfsrcread[0].oldsurf, macrobodies)

    def __str__(self):
        card = 'EXTERNALSOURCE\n'
        for source in self._source:
            card += str(source)
        for surfsrcread in self._surfsrcread:
            card += str(surfsrcread)
        for distribution in self._distributions:
            card += str(distribution)
        card += self._unparsed
        card += '\n\n'
        return card


class Source(BaseModel):

    ###########有问题的card_option_types################################
    # card_option_types = {
    #     'SOURCE': [int],
    #     'FRACTION': [float],
    #     'PARTICLE': [int],
    #     'POINT': ['list', float],
    #     'SPHERE': ['list', str],
    #     'CYL/X': ['list', str],
    #     'CYL/Y': ['list', str],
    #     'CYL/Z': ['list', str],
    #     'SURFACE': ['list', str],
    #     'X': [str],
    #     'Y': [str],
    #     'Z': [str],
    #     'POSITION': ['list', str],
    #     'RADIUS': [str],
    #     'AXIS': ['list', str],
    #     'POlAR': ['list', str],
    #     'POLARTHETA': ['list', str],
    #     'EXTENT': ['list', str],
    #     'HEIGHT': ['list', str],
    #     'NORM': ['list', str],
    #     'VECTOR': ['list', str],
    #     'DIRECMIU': ['list', str],
    #     'FAIVECTOR': ['list', str],
    #     'DIRECFAI': ['list', str],
    #     'ENERGY': ['list', str],
    #     'CELL': ['list', str],
    #     'SAMPEFF': [str],
    #     'BIASFRAC': [str],
    #     'WEIGHT': [str],
    #     'TRANSFORM': [str]
    # }

    def __init__(self, source_id=None, fraction=None, particle=None, point=None, sphere=None, cyl_x=None, cyl_y=None,
                 cyl_z=None, surface=None, x=None, y=None, z=None, position=None, radius=None, axis=None, polar=None,
                 polartheta=None, extent=None, height=None, norm=None, vector=None, direcmiu=None, faivector=None,
                 direcfai=None, energy=None, cell=None, sampeff=None, biasfrac=None, weight=None, transform=None):
        self._id = source_id
        self._fraction = fraction
        self._particle = particle
        self._point = point
        self._sphere = sphere
        self._cyl_x = cyl_x
        self._cyl_y = cyl_y
        self._cyl_z = cyl_z
        self._surface = surface
        self._x = x
        self._y = y
        self._z = z
        self._position = position
        self._radius = radius
        self._axis = axis
        self._polar = polar
        self._polartheta = polartheta
        self._extent = extent
        self._height = height
        self._norm = norm
        self._vector = vector
        self._direcmiu = direcmiu
        self._faivector = faivector
        self._direcfai = direcfai
        self._energy = energy
        self._cell = cell
        self._sampeff = sampeff
        self._biasfrac = biasfrac
        self._weight = weight
        self._transform = transform

    def add_options(self, options):
        if 'SOURCE' in options.keys():
            self._id = options['SOURCE']
        if 'FRACTION' in options.keys():
            self._fraction = options['FRACTION']
        if 'PARTICLE' in options.keys():
            self._particle = options['PARTICLE']
        if 'POINT' in options.keys():
            self._point = options['POINT']
        if 'SPHERE' in options.keys():
            self._sphere = options['SPHERE']
        if 'CYL/X' in options.keys():
            self._cyl_x = options['CYL/X']
        if 'CYL/Y' in options.keys():
            self._cyl_y = options['CYL/Y']
        if 'CYL/Z' in options.keys():
            self._cyl_z = options['CYL/Z']
        if 'SURFACE' in options.keys():
            self._surface = options['SURFACE']
        if 'X' in options.keys():
            self._x = options['X']
        if 'Y' in options.keys():
            self._y = options['Y']
        if 'Z' in options.keys():
            self._z = options['Z']
        if 'POSITION' in options.keys():
            self._position = options['POSITION']
        if 'RADIUS' in options.keys():
            self._radius = options['RADIUS']
        if 'AXIS' in options.keys():
            self._axis = options['AXIS']
        if 'EXTENT' in options.keys():
            self._extent = options['EXTENT']
        if 'POLAR' in options.keys():
            self._polar = options['POLAR']
        if 'POLARTHETA' in options.keys():
            self._polartheta = options['POLARTHETA']
        if 'HEIGHT' in options.keys():
            self._height = options['HEIGHT']
        if 'NORM' in options.keys():
            self._norm = options['NORM']
        if 'VECTOR' in options.keys():
            self._vector = options['VECTOR']
        if 'DIRECMIU' in options.keys():
            self._direcmiu = options['DIRECMIU']
        if 'FAIVECTOR' in options.keys():
            self._faivector = options['FAIVECTOR']
        if 'DIRECFAI' in options.keys():
            self._direcfai = options['DIRECFAI']
        if 'ENERGY' in options.keys():
            self._energy = options['ENERGY']
        if 'CELL' in options.keys():
            self._cell = options['CELL']
        if 'SAMPEFF' in options.keys():
            self._sampeff = options['SAMPEFF']
        if 'BIASFRAC' in options.keys():
            self._biasfrac = options['BIASFRAC']
        if 'WEIRGT' in options.keys():
            self._weight = options['WEIRGT']
        if 'TRANSFORM' in options.keys():
            self._transform = options['TRANSFORM']

    def check(self):
        assert self._id > 0

    def postprocess(self, geometry=None):
        """
        处理通用源中不采用Distribution的Cell部分，使其兼容MCNP风格的栅元展开式
        """
        if self._cell is None:
            return

        value_str = ' '.join([str(x) for x in self._cell])
        if not re.search(r'\(|\)', value_str):
            return

        cell_expansion = ' '.join([str(x) for x in self._cell])
        post_value = []
        Distribution.process_cell_expansion(expansion=cell_expansion, geometry=geometry, value=post_value)
        self._cell = post_value

    @property
    def surface(self):
        return self._surface

    # def __str__(self):
    #     card = 'Source'
    #     card += ' ' + str(self._id)
    #     if self._fraction is not None:
    #         card += ' Fraction = ' + str(self._fraction)
    #     # if self._particle is not None:
    #     #     card += ' Particle = ' + ' '.join([str(x) for x in self._particle])
    #     if self._particle is not None:
    #         card += ' Particle = ' + str(self._particle)
    #     if self._point is not None:
    #         card += ' point = ' + ' '.join([str(x) for x in self._point])
    #     if self._sphere is not None:
    #         card += ' sphere = ' + ' '.join([str(x) for x in self._sphere])
    #     if self._cyl_x is not None:
    #         card += ' Cyl/x = ' + ' '.join([str(x) for x in self._cyl_x])
    #     if self._cyl_y is not None:
    #         card += ' Cyl/y = ' + ' '.join([str(x) for x in self._cyl_y])
    #     if self._cyl_z is not None:
    #         card += ' Cyl/z = ' + ' '.join([str(x) for x in self._cyl_z])
    #     if self._surface is not None:
    #         card += ' Surface = ' + ' '.join([str(x) for x in self._surface])
    #     # if self._x is not None:
    #     #     card += ' X = ' + ' '.join([str(x) for x in self._x])
    #     # if self._y is not None:
    #     #     card += ' Y = ' + ' '.join([str(x) for x in self._y])
    #     # if self._z is not None:
    #     #     card += ' Z = ' + ' '.join([str(x) for x in self._z])
    #     if self._x is not None:
    #         card += ' X = ' + str(self._x)
    #     if self._y is not None:
    #         card += ' Y = ' + str(self._y)
    #     if self._z is not None:
    #         card += ' Z = ' + str(self._z)
    #     if self._position is not None:
    #         card += ' Position = ' + ' '.join([str(x) for x in self._position])
    #     # if self._radius is not None:
    #     #     card += ' Radius = ' + ' '.join([str(x) for x in self._radius])
    #     if self._radius is not None:
    #         card += ' Radius = ' + str(self._radius)
    #     if self._axis is not None:
    #         card += ' Axis = ' + ' '.join([str(x) for x in self._axis])
    #     if self._extent is not None:
    #         card += ' Extent = ' + ' '.join([str(x) for x in self._extent])
    #     if self._polar is not None:
    #         card += ' Polar = ' + ' '.join([str(x) for x in self._polar])
    #     if self._polartheta is not None:
    #         card += ' PolarTheta = ' + ' '.join([str(x) for x in self._polartheta])
    #     if self._height is not None:
    #         card += ' Height = ' + ' '.join([str(x) for x in self._height])
    #     if self._norm is not None:
    #         card += ' Norm = ' + ' '.join([str(x) for x in self._norm])
    #     if self._vector is not None:
    #         card += ' Vector = ' + ' '.join([str(x) for x in self._vector])
    #     if self._direcmiu is not None:
    #         card += 'DirecMiu = ' + ' '.join([str(x) for x in self._direcmiu])
    #     if self._faivector is not None:
    #         card += ' FaiVector = ' + ' '.join([str(x) for x in self._faivector])
    #     if self._direcfai is not None:
    #         card += ' DirecFai = ' + ' '.join([str(x) for x in self._direcfai])
    #     if self._energy is not None:
    #         card += ' Energy = ' + ' '.join([str(x) for x in self._energy])
    #     if self._cell is not None:
    #         card += ' Cell = ' + ' '.join([str(x) for x in self._cell])
    #     if self._sampeff is not None:
    #         card += ' SampEff = ' + str(self._sampeff)
    #     if self._biasfrac is not None:
    #         card += ' Biasfrac = ' + str(self._biasfrac)
    #     if self._weight is not None:
    #         card += ' Weight = ' + str(self._weight)
    #     if self._transform is not None:
    #         card += ' Transform= ' + str(self._transform)
    #     card += '\n'
    #     return card
    def __str__(self):
        card = 'Source ' + str(self._id)
        if self._fraction is not None:
            card += ' Fraction=' + str(self._fraction)
        
        # 对于每个属性，检查是否是列表类型，然后相应地处理
        attributes = [
            ('Particle', self._particle),
            ('point', self._point),
            ('sphere', self._sphere),
            ('Cyl/x', self._cyl_x),
            ('Cyl/y', self._cyl_y),
            ('Cyl/z', self._cyl_z),
            ('Surface', self._surface),
            ('X', self._x),
            ('Y', self._y),
            ('Z', self._z),
            ('Position', self._position),
            ('Radius', self._radius),
            ('Axis', self._axis),
            ('Extent', self._extent),
            ('Polar', self._polar),
            ('PolarTheta', self._polartheta),
            ('Height', self._height),
            ('Norm', self._norm),
            ('Vector', self._vector),
            ('DirecMiu', self._direcmiu),
            ('FaiVector', self._faivector),
            ('DirecFai', self._direcfai),
            ('Energy', self._energy),
            ('Cell', self._cell),
            ('SampEff', self._sampeff),
            ('Biasfrac', self._biasfrac),
            ('Weight', self._weight),
            ('Transform', self._transform),
        ]
        
        for name, value in attributes:
            if value is not None:
                # 如果是列表，则使用 join，否则直接转换为字符串
                card_part = ' '.join(map(str, value)) if isinstance(value, list) else str(value)
                card += ' ' + name + '=' + card_part
        
        card += '\n'
        return card


class Distribution(BaseModel):
    card_option_types = {
        'DEPEND': [str],
        'TYPE': [int],
        'VALUE': ['list', float, -1],
        'PROBABILITY': ['list', float, -1],
        'BIAS': ['list', float, -1]
    }

    def __init__(self, id=None, depend=None, type=None, value=None, probability=None, bias=None):
        self._id = id
        self._depend = depend
        self._type = type
        self._value = value
        self._probability = probability
        self._bias = bias

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, _id):
        self._id = _id

    @property
    def depend(self):
        return self._depend

    @depend.setter
    def depend(self, d):
        self._depend = d

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, t):
        self._type = t

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, v):
        self._value = v

    @property
    def probability(self):
        return self._probability

    @probability.setter
    def probability(self, prob):
        self._probability = prob

    @property
    def bias(self):
        return self._bias

    @bias.setter
    def bias(self, b):
        self._bias = b

    @staticmethod
    def is_integral(number):
        return isinstance(number, float) and int(number) == number

    @staticmethod
    def process_cell_expansion(expansion=None, geometry=None, value=None):
        """
        将MCNP风格的栅元展开式，修改为RMC风格的栅元展开式
        """
        if not re.search(r'\(|\)', expansion):
            return
        # split the single expansion, ['2', ' 2  1 1', ' 2 3 4 ', ' 6']
        cell_expansion = [x for x in re.split(r'>|\(|\)', expansion) if x != ' ']
        scope = []
        fill = []
        while cell_expansion:
            cell_lattice = re.findall(r'\d+', cell_expansion[0])

            if len(cell_lattice) == 1:
                cell = geometry.cell_dict[int(cell_lattice[0])]
                # for the bottom cell without filling
                if cell.fill is None:
                    value.append(int(cell_lattice[0]))
                    del cell_expansion[0]
                    continue
                else:
                    # for cell with filling
                    fill_universe = geometry.univ_dict[cell.fill]
                    # for universe with lattice, get lattice scope
                    if fill_universe.lattice is not None:
                        scope = fill_universe.lattice.scope
                        fill = fill_universe.lattice.fill

                    value.append(int(cell_lattice[0]))
                    value.append('>')
                    del cell_expansion[0]

            elif len(cell_lattice) == 3:
                x = int(cell_lattice[0])
                y = int(cell_lattice[1])
                z = int(cell_lattice[2])

                lattice_index = x + scope[0] * (y - 1) + scope[0] * scope[1] * (z - 1)

                universe_id = fill[lattice_index - 1]
                fill_universe = geometry.univ_dict[universe_id]
                if fill_universe.lattice is not None:
                    scope = fill_universe.lattice.scope
                    fill = fill_universe.lattice.fill

                value.append(lattice_index)
                value.append('>')
                del cell_expansion[0]

    def postprocess(self, geometry=None):
        """
        兼容MCNP风格的通用源Cell输入，即采用三点坐标表示lattice位置，如2 > (2 1 1)> > 6，并将其改为
        RMC风格的通用源输入，即使用数字序号表示lattice文职，如2 > 3 > 6 (2*2 lattice)
        """
        if not self.value:
            return
        value_str = ' '.join([str(int(x)) if Distribution.is_integral(x) else str(x) for x in self.value])

        if not re.search(r'>|\(|\)', value_str):
            return
        # get all the cell expansion with (), for example: 
        # expansion: 21 > 5>6>60 2 > ( 2  1 1)>( 2 3 4 ) >  6  2 > (1 1 1)>( 2 3 5)> 66 
        # get [2 > ( 2  1 1)>( 2 3 4 ) >  6, 2 > (1 1 1)>( 2 3 5)> 66]
        cell_compile = re.compile(r'\d+\s*(?:>\s*(?:\d+|\([\d\s]+\)\s*))*')
        cell_list = re.findall(cell_compile, value_str)
        # replace the expansion using the lattice index
        post_value = []
        for expansion in cell_list:
            self.process_cell_expansion(expansion=expansion, geometry=geometry, value=post_value)
        self._value = post_value

    def __str__(self):
        card = "Distribution "
        card += str(self.id)
        if self.depend:
            card += ' depend = ' + str(self.depend)
        card += ' type = ' + str(self.type)
        if self.value:
            card += ' value = ' + ' '.join(
                [str(int(x)) if Distribution.is_integral(x) else str(x) for x in self.value])
        if self.probability:
            card += ' probability = ' + ' '.join([f'{x:.12g}' for x in self.probability])
        if self.bias:
            card += ' bias = ' + ' '.join([f'{x:.12g}' for x in self.bias])
        card += '\n'
        return card


class Surfsrcread(BaseModel):
    card_option_types = {
        'OLDSURF': ['list', float],
        'NEWSURF': ['list', float],
        'CELL': ['list', int],
        'PARTYPE': [int],
        'COLI': [int],
        'WTM': [float],
        'AXIS': ['list', float],
        'EXTENT': ['list', float],
        'POSACE': ['list', float],
        'TR': ['list', int]
    }

    def __init__(self, oldsurf=None, newsurf=None, cell=None, partype=None, coli=None, wtm=None, axis=None, extent=None,
                 posace=None, tr=None):
        self._oldsurf = oldsurf
        self._newsurf = newsurf
        self._cell = cell
        self._partype = partype
        self._coli = coli
        self._wtm = wtm
        self._axis = axis
        self._extent = extent
        self._posace = posace
        self._tr = tr

    def check(self):
        pass

    @property
    def newsurf(self):
        return self._newsurf

    @property
    def oldsurf(self):
        return self._oldsurf

    def add_options(self, options):
        if 'OLDSURF' in options.keys():
            self._oldsurf = options['OLDSURF']
        if 'NEWSURF' in options.keys():
            self._newsurf = options['NEWSURF']
        if 'CELL' in options.keys():
            self._cell = options['CELL']
        if 'PARTYPE' in options.keys():
            self._partype = options['PARTYPE']
        if 'COLI' in options.keys():
            self._coli = options['COLI']
        if 'WTM' in options.keys():
            self._wtm = options['WTM']
        if 'AXIS' in options.keys():
            self._axis = options['AXIS']
        if 'EXTENT' in options.keys():
            self._extent = options['EXTENT']
        if 'POSACE' in options.keys():
            self._posace = options['POSACE']
        if 'TR' in options.keys():
            self._tr = options['TR']

    def __str__(self):
        card = 'SurfSrcRead'
        if self._oldsurf is not None:
            card += ' OldSurf =' + ' '.join([str(x) for x in self._oldsurf])
        if self._newsurf is not None:
            card += ' NewSurf = ' + ' '.join([str(x) for x in self._newsurf])
        if self._cell is not None:
            card += ' Cell = ' + ' '.join([str(x) for x in self._cell])
        if self._partype is not None:
            card += ' Partype = ' + str(self._partype)
        if self._coli is not None:
            card += ' Coli = ' + str(self._coli)
        if self._wtm is not None:
            card += ' Wtm = ' + str(self._wtm)
        if self._axis is not None:
            card += ' Axis = ' + ' '.join([str(x) for x in self._axis])
        if self._extent is not None:
            card += ' Extent = ' + ' '.join([str(x) for x in self._extent])
        if self._posace is not None:
            card += ' Posace = ' + ' '.join([str(x) for x in self._posace])
        if self._tr is not None:
            card += ' Tr = ' + ' '.join([str(x) for x in self._tr])
        card += '\n'
        return card
