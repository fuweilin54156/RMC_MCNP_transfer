# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2019-11-15


from RMC.model.input.base import YMLModelObject as BaseModel
import numpy as np
import warnings
import copy
import re
from RMC.model.input.Surface import *


class Cell(BaseModel):
    yaml_tag = u'!cell'
    card_option_types = {
        'FILL': [int],
        'CELL': [int],
        'MAT': ['list', int, -1],
        'VOL': [float],
        'TMP': [float],
        'DENS': [float],
        'VOID': [bool],
        'INNER': [bool],
        'MOVE': ['list', float, 3],
        'ROTATE': ['list', float, 9],
        'NOBURN': [bool],
        'IMP:N': [int],
        'IMP:P': [int],
        'IMP:E': [int]
    }

    def __init__(self, name=None, number=-1, bounds='', material=None, volume=None, fill=None, inner=False,
                 temperature=None, density=None, void=False, transformation=None, noburn=0, imp_n=None,
                 imp_p=None, imp_e=None):
        self.name = name
        self.number = number
        self.bounds = bounds
        self.fill = fill
        self.inner = inner
        self.material = material
        self.temperature = temperature
        self.volume = volume
        self.density = density
        self.void = void
        self.transformation = transformation
        self.include = None  # included universe
        self.noburn = noburn
        self.imp_n = imp_n
        self.imp_p = imp_p
        self.imp_e = imp_e

    def check(self):
        # assert self.temperature > 0
        assert self.number >= 0

    def add_bounds(self, bounds):
        self.bounds = bounds
        pass

    def add_options(self, options):
        self.fill = options['FILL']
        self.material = options['MAT']
        self.volume = options['VOL']
        self.temperature = options['TMP']
        self.density = options['DENS']
        self.void = options['VOID']
        self.inner = options['INNER']
        self.noburn = options['NOBURN']
        self.imp_n = options['IMP:N']
        self.imp_p = options['IMP:P']
        self.imp_e = options['IMP:E']
        if options['MOVE'] is not None or options['ROTATE'] is not None:
            self.transformation = Transformation(move=options['MOVE'], rotate=options['ROTATE'])

    def count_cell(self, cell_id):
        if self.number == cell_id:
            return 1
        if self.fill is None:
            return 0
        elif self.include is None:
            raise ValueError("Postprocessing for cell {} has not been done!".format(self.number))
        else:
            return self.include.count_cell(cell_id)

    def duplicate(self, geom, end_univ=None, excepts: list = None, template_univ=None, default_pin=None):
        new_end_univ = None
        new_cell = Cell(name=self.name, number=geom.get_new_cell_id(), bounds=self.bounds, material=self.material,
                        volume=self.volume, fill=self.fill, inner=self.inner, temperature=self.temperature,
                        density=self.density, void=self.void, transformation=copy.copy(self.transformation))
        geom.add_cell(new_cell)
        if new_cell.fill is not None:
            new_univ, new_end_univ = geom.get_univ(new_cell.fill).duplicate(geom=geom, end_univ=end_univ,
                                                                            excepts=excepts,
                                                                            template_univ=template_univ,
                                                                            default_pin=default_pin)
            new_cell.fill = new_univ.number
        return new_cell, new_end_univ

    def __str__(self):
        s = f'CELL {self.number} {self.bounds} '
        if self.material is not None:
            s += 'mat=' + ''.join([str(x) for x in self.material])
            s += ' '
        if self.fill is not None:
            s += f'fill={self.fill} '
        if self.inner:
            s += 'inner=1 '
        if self.temperature is not None:
            s += f'tmp={self.temperature} '
        if self.volume is not None:
            s += f'vol={self.volume} '
        if self.density is not None:
            s += f'dens={self.density} '
        if self.void:
            s += 'void=1 '
        if self.noburn:
            s += 'noburn=1'
        if self.imp_n:
            s += f'imp:n={self.imp_n} '
        if self.imp_p:
            s += f'imp:p={self.imp_p} '
        if self.imp_e:
            s += f'imp:e={self.imp_e} '
        if self.transformation is not None:
            s += str(self.transformation) + ' '
        return s.strip() + '\n'


class Universe(BaseModel):
    yaml_tag = u'!universe'
    card_option_types = {
        'FILL': ['list', int, -1],
        'SCOPE': ['list', int, 3],
        'PITCH': ['list', float, 3],
        'LAT': [int],
        'MOVE': ['list', float, 3],
        'ROTATE': ['list', float, 9],
        'SITA': [float],
        # only for lattice type 5
        'MAX': [int],
        'LENGTH': [float],
        'RADIUS': [float],
        'BALLUNI': [int],
        'INTERVALUNI': [int],
        # for random geometry
        'MATRIC': [int],
        'PARTICLE': ['list', int, -1],
        'PF': ['list', float, -1],
        'RAD': ['list', float, -1],
        'RSA': [int],
        'TYPE': [int],
        'SIZE': ['list', float, -1],
        'DEM': [int],
        'TIME': [float],
        'PFCORRECT': [int]
    }

    def __init__(self, name=None, number=-1, cells=None, transformation=None, lattice=None):
        if cells is None:
            cells = []
        self.name = name
        self.number = number
        self.cells = cells
        self.transformation = transformation
        self.lattice = lattice
        self.include = set()

    def check(self):
        assert self.number >= 0

    def add_options(self, options):
        if options['MOVE'] is not None or options['ROTATE'] is not None:
            self.transformation = Transformation(move=options['MOVE'], rotate=options['ROTATE'])
        if options['LAT'] is not None:
            self.lattice = Lattice(type=options['LAT'], pitch=options['PITCH'], scope=options['SCOPE'],
                                   fill=options['FILL'], theta=options['SITA'], max_scope=options['MAX'],
                                   fill_ball=options['BALLUNI'], fill_interval=options['INTERVALUNI'],
                                   radius=options['RADIUS'], length=options['LENGTH'], matric=options['MATRIC'],
                                   particle=options['PARTICLE'], pf=options['PF'], rad=options['RAD'],
                                   rsa=options['RSA'], typerg=options['TYPE'], sizerg=options['SIZE'],
                                   dem=options['DEM'], time=options['TIME'], pfcorrect=options['PFCORRECT'])

    def add_cells(self, cells):
        self.cells.extend(cells)

    def add_cell(self, cell):
        self.cells.append(cell)

    def postprocess(self):
        # todo: decorate the postprocess method with checking.
        self.check()
        for cell in self.cells:
            cell.postprocess()
        if self.lattice is not None:
            self.lattice.postprocess()
        if self.transformation is not None:
            self.transformation.postprocess()

    def count_cell(self, cell_id):
        num = 0
        if self.lattice is None:
            for cell in self.cells:
                num += cell.count_cell(cell_id)
        elif self.include is None:
            raise ValueError("Postprocessing for universe {} has not been done!".format(self.number))
        else:
            # todo: refactor needed.
            unique, counts = np.unique(self.lattice.fill, return_counts=True)
            occurrence = dict(zip(unique, counts))
            for univ in self.include:
                if univ.number in occurrence:
                    num += occurrence[univ.number] * univ.count_cell(cell_id)
        return num

    # todo: only 2D positions in the excepts list is allowed currently.
    # todo: hash code of the lattice list can be used to accelerate.
    def compare_pattern(self, univ, excepts: list = None, default_pin=None):
        if self.lattice is None or univ.lattice is None:
            raise ValueError("The two universes do not both have lattices.")
        return self.lattice.compare_pattern(univ.lattice, excepts=excepts, default_pin=default_pin)

    def duplicate(self, geom, end_univ=None, excepts: list = None, template_univ=None, default_pin=None):
        new_univ = Universe(name=self.name, number=geom.get_new_univ_id(), cells=[],
                            transformation=copy.copy(self.transformation), lattice=copy.copy(self.lattice))
        new_end_univ = None
        geom.add_universe(new_univ, postprocessing=True)

        # todo: end_univ and excepts.
        if end_univ is not None and self.number == end_univ.number:
            if excepts is None or excepts == []:
                return new_univ, new_end_univ
            elif self.lattice is None:
                raise TypeError("Duplicates with excepts should end with a universe with lattice.")
            elif template_univ is None and default_pin is None:
                raise TypeError("Template universe or default pin should be specified.")
            elif template_univ is not None and default_pin is not None:
                raise TypeError("Template universe and default pin should be both specified.")
            elif template_univ is not None and template_univ.lattice is None:
                raise TypeError("The template universe should be with lattice.")
            else:
                new_end_univ = new_univ
                if template_univ is not None:
                    new_univ.lattice.sub(excepts=excepts, template_lattice=template_univ.lattice)
                else:
                    new_univ.lattice.sub(excepts=excepts, default_pin=default_pin)

        for idx in range(len(self.cells)):
            new_cell, new_end_univ_temp = self.cells[idx].duplicate(geom=geom, end_univ=end_univ, excepts=excepts,
                                                                    template_univ=template_univ,
                                                                    default_pin=default_pin)
            new_univ.add_cell(new_cell)
            if new_end_univ_temp is not None:
                new_end_univ = new_end_univ_temp
        if self.lattice is not None and (end_univ is None or self.number != end_univ.number):
            mapping = {}
            for univ in sorted(self.include, key=lambda x: x.number):
                new_sub_univ, new_end_univ_temp = univ.duplicate(geom=geom, end_univ=end_univ, excepts=excepts,
                                                                 template_univ=template_univ, default_pin=default_pin)
                if new_end_univ_temp is not None:
                    new_end_univ = new_end_univ_temp
                mapping[univ.number] = new_sub_univ.number
            for idx, val in enumerate(new_univ.lattice.fill):
                new_univ.lattice.fill[idx] = mapping[val]
        return new_univ, new_end_univ

    # todo: refactor the codes.
    def find_next_lattice_univ(self):
        if self.lattice is not None:
            return self.number
        for cell in self.cells:
            if cell.fill is not None:
                univ = cell.include.find_next_lattice_univ()
                if univ is not None:
                    return univ
        return None

    def find_contain_univ(self, univ_set: set):
        if self.number in univ_set:
            return self.number
        if self.lattice is None:
            for cell in self.cells:
                if cell.fill is not None:
                    uid = cell.include.find_contain_univ(univ_set)
                    if uid is not None:
                        return uid
            return None
        else:
            for univ in self.include:
                uid = univ.find_contain_univ(univ_set)
                if uid is not None:
                    return uid
            return None

    def __str__(self):
        s = f'UNIVERSE {self.number} '
        if self.transformation is not None:
            s += str(self.transformation) + ' '
        if self.lattice is not None:
            s += str(self.lattice)
        s += '\n'
        for cell in self.cells:
            s += str(cell)
        s += '\n'
        return s

    def __iter__(self):
        for cell in self.cells:
            yield cell


class Lattice(BaseModel):
    yaml_tag = u'!lattice'

    def __init__(self, pitch=np.zeros(3), scope=np.zeros(3), fill=np.ones(1), type=None, theta=None,
                 length=None, radius=None, fill_ball=None, fill_interval=None, max_scope=None,  # lat=5
                 matric=None, particle=None, pf=None, rad=None, typerg=None, sizerg=None, pfcorrect=False, rsa=None,
                 dem=None, time=None  # lat=3 4
                 ):
        self.pitch = np.array(pitch)
        self.scope = np.array(scope)
        self.fill = np.array(fill)
        self.type = type
        self.theta = theta
        # lat = 3 4
        self.matric = matric
        self.particle = particle
        self.pf = pf
        self.rad = rad
        self.typerg = typerg
        self.sizerg = sizerg
        self.pfcorrect = pfcorrect
        self.rsa = rsa
        self.dem = dem
        self.time = time
        # lat = 5
        self.length = length
        self.radius = radius
        self.fill_ball = fill_ball
        self.fill_interval = fill_interval
        self.max_scope = max_scope

    def check(self):
        assert self.type >= 0

    def __str__(self):
        if self.type in [1, 2]:
            s = f'LAT={self.type} SCOPE='
            for ele in self.scope:
                s += f'{ele} '
            s += 'PITCH='
            for ele in self.pitch:
                s += f'{ele} '
            if self.theta is not None and self.type == 2:
                s += f'SITA={self.theta} '
            s += 'FILL='
            if len(set(self.fill)) == 1:
                s += f'{self.fill[0]} '
                s += f'* {len(self.fill)} \n'
            else:
                for idx, val in enumerate(self.fill):
                    if idx % self.scope[0] == 0:
                        s += '\n '
                    s += f'{val} '
            return s.replace(' \n', '\n').strip()
        elif self.type in [3, 4]:
            s = f'LAT={self.type} '

            s += f'MATRIC={self.matric} '
            s += 'PARTICLE='
            for ele in self.particle:
                s += f'{ele} '
            s += 'PF='
            for ele in self.pf:
                s += f'{ele} '
            s += 'RAD='
            for ele in self.rad:
                s += f'{ele} '
            s += ' TYPE=%d ' % self.typerg
            s += 'SIZE='
            for ele in self.sizerg:
                s += f'{ele} '
            if self.pfcorrect is not None:
                s += f' PFCORRECT={self.pfcorrect} '
            if self.rsa is not None:
                s += f' RSA={self.rsa} '
            if self.dem is not None:
                s += f' DEM={self.dem} '
            if self.time is not None:
                s += f' TIME={self.time} '
        elif self.type == 5:
            s = f'LAT={self.type} \n'
        return s.replace(' \n', '\n').strip()

    def shape_equal(self, other):
        return np.all(self.pitch == other.pitch) and np.all(self.scope == other.scope) \
               and self.fill.shape == other.fill.shape and self.type == self.type and self.theta == other.theta

    def compare_pattern(self, other, excepts=None, default_pin=None):
        if not self.shape_equal(other):
            raise ValueError("The two lattices have different sizes.")
        except_pos = {}
        if excepts is not None:
            except_idx = 0
            for pair in excepts:
                except_pos[pair[0] * self.scope[1] + pair[1]] = except_idx
                except_idx += 1

        pattern_same = True
        except_same = True
        for idx, val in enumerate(other.fill):
            if idx in except_pos:
                compare_val = self.fill[idx]
                if default_pin is not None:
                    compare_val = default_pin[except_pos[idx]]
                if val != compare_val:
                    except_same = False
            else:
                if val != self.fill[idx]:
                    pattern_same = False
            if not except_same and not pattern_same:
                break

        return [pattern_same, except_same]

    def sub(self, excepts, template_lattice=None, default_pin=None):
        if (template_lattice is None and default_pin is None) or \
           (template_lattice is not None and default_pin is not None):
            raise TypeError("One and only one of the template lattice and default pin should be specified.")
        idx = 0
        for pair in excepts:
            pos = pair[0] * self.scope[1] + pair[1]
            if template_lattice is not None:
                self.fill[pos] = template_lattice.fill[pos]
            else:
                if default_pin[idx] is not None:
                    self.fill[pos] = default_pin[idx]
            idx += 1

    def get_idx(self, position):
        if len(position) == 2:
            return position[0] * self.scope[1] + position[1]
        elif len(position) == 3:
            return (position[0] * self.scope[1] + position[1]) * self.scope[2] + position[2]
        else:
            raise ValueError("Lattice position should be 2-dim or 3-dim.")

    def __copy__(self):
        return Lattice(pitch=np.copy(self.pitch), scope=np.copy(self.scope), fill=np.copy(self.fill), type=self.type,
                       theta=self.theta)


class Transformation(BaseModel):
    yaml_tag = u'!transformation'

    def __init__(self, move=None, rotate=None):
        if move is not None:
            self.move = np.array(move)
        else:
            self.move = None
        if rotate is not None:
            self.rotate = np.array(rotate)
        else:
            self.rotate = None

    def check(self):
        assert self.move.shape == tuple([3]) and self.rotate.shape == tuple([3, 3])

    def __str__(self):
        s = ''
        if self.move is not None:
            s += 'move='
            for idx, val in enumerate(self.move):
                s += f'{val:.12g} '
        if self.rotate is not None:
            s += 'rotate='
            for idx, val in enumerate(self.rotate):
                s += '{:.12g} '.format(val)
        return s.strip()

    def __copy__(self):
        move = None
        if self.move is not None:
            move = np.copy(self.move)
        rotate = None
        if self.rotate is not None:
            rotate = np.copy(self.rotate)
        return Transformation(move=move, rotate=rotate)


class Geometry(BaseModel):
    yaml_tag = u'!geometry'

    def __init__(self, universes=None):
        if universes is None:
            universes = []
        self.universes = universes
        self.univ_dict = {}
        self.cell_dict = {}

    def check(self):
        assert len(self.universes) > 0

    def add_universe(self, univ, postprocessing=False):
        self.universes.append(univ)
        if postprocessing:
            self.univ_dict[univ.number] = univ
            self.postprocess_for(univ)

    def add_cell(self, cell, postprocesing=True):
        if postprocesing:
            self.cell_dict[cell.number] = cell
            if cell.fill is not None:
                cell.include = self.get_univ(cell.fill)
                self.postprocess_for(self.get_univ(cell.fill))

    def postprocess(self):
        for univ in self.universes:
            self.univ_dict[univ.number] = univ
            for cell in univ.cells:
                self.cell_dict[cell.number] = cell
        for univ in self.universes:
            univ.postprocess()
            if univ.lattice is not None and univ.lattice.type in [1, 2]:
                for uid in univ.lattice.fill:
                    univ.include.add(self.univ_dict[uid])
            for cell in univ.cells:
                if cell.fill is not None:
                    cell.include = self.univ_dict[cell.fill]

    def postprocess_for(self, univ=None):
        if univ is None:
            self.postprocess()
            return
        univ.postprocess()
        if univ.lattice is not None:
            for uid in univ.lattice.fill:
                univ.include.add(self.univ_dict[uid])
        # todo: double postprocess exists.
        for cell in univ.cells:
            if cell.fill is not None:
                cell.include = self.univ_dict[cell.fill]

    # 对栅元的布尔运算修改
    def trans_cell(self, macrobodies=None):
        for univ in self.universes:
            for cell in univ.cells:
                # bounds = re.compile(r'[^&:]?-?!?\d+\.?\d*').findall(cell.bounds)
                # cell.bounds = re.compile(r'[&]+|[:]+|[()]|-?!?\d+\.?\d*').findall(cell.bounds)  # 匹配栅元中的布尔表达式
                # 匹配栅元中的布尔表达式,提取包含!+数字和所有数字[(!1,2:,3)&,-1.3&,4.2:]这种list
                cell.bounds = re.compile(r'[\(\s]*-?!?\d+\.?\d*[\s\)]*\s?&?\s?:?').findall(cell.bounds)
                for i, card in enumerate(cell.bounds):
                    if '!' in card:
                        continue
                    else:
                        bound = re.findall('-?\d+\.?\d*', card)[0]  # 匹配得到每个字符串的数字
                        if abs(float(bound)) - int(abs(float(bound))) == 0:
                            if abs(float(bound)) not in macrobodies.keys():
                                continue
                            else:
                                body = macrobodies[abs(float(bound))]
                                cell.bounds[i] = re.sub(bound, body.bound(in_=float(bound) < 0), card)
                        else:
                            body_id = abs(int(bound.split('.')[0]))
                            surface_id = int(bound.split('.')[1])
                            body = macrobodies[body_id]
                            surface_id = body.surf_ids[surface_id - 1]
                            if float(bound) < 0:
                                cell.bounds[i] = re.sub(bound, str(-surface_id), card)
                            else:
                                cell.bounds[i] = re.sub(bound, str(surface_id), card)
                cell.bounds = ''.join(cell.bounds)

    def get_univ(self, uid):
        return self.univ_dict[uid]

    def get_cell(self, cid):
        return self.cell_dict[cid]

    def proc_lat5(self, surf_id):
        surf_max_id = int(surf_id)
        surfs = []
        for univ in self.universes:
            if univ.lattice:
                if univ.lattice.type == 5:
                    cell_max_id = self.get_max_cell_id()
                    uni_max_id = self.get_max_uni_id()
                    # get the processed cells, surfs and univ id
                    cell_num_list = [cell_max_id + i for i in range(1, 16)]
                    univ_num = uni_max_id + 1
                    surf_num_list = [surf_max_id + i for i in range(1, 15)]
                    surf_max_id += 15  # there are around 15 additional cells in processing one lat 5 structure

                    # Lat process
                    univ.lattice.length = univ.lattice.length * np.sqrt(2)
                    post_scope = np.array([univ.lattice.max_scope, univ.lattice.max_scope, univ.lattice.max_scope])
                    post_pitch = np.array([univ.lattice.length, univ.lattice.length, univ.lattice.length])
                    post_fill = np.ones(univ.lattice.max_scope * univ.lattice.max_scope * univ.lattice.max_scope,
                                        dtype=int) * univ_num
                    post_lat = Lattice(type=1, scope=post_scope, pitch=post_pitch, fill=post_fill)

                    # Univ process
                    post_univ = Universe(number=univ_num)  # the processed Universe for Lat = 5

                    # Surfs process
                    l = univ.lattice.length
                    r = univ.lattice.radius
                    param_list = [[0, 0, 0, r], [l, 0, 0, r], [l, l, 0, r], [0, l, 0, r], [0, 0, l, r], [l, 0, l, r],
                                  [l, l, l, r], [0, l, l, r], [l / 2, 0, l / 2, r], [l, l / 2, l / 2, r],
                                  [l / 2, l, l / 2, r], [0, l / 2, l / 2, r], [l / 2, l / 2, l, r],
                                  [l / 2, l / 2, 0, r]]
                    for i in range(14):
                        surfs.append(Surface(number=surf_num_list[i], parameters=param_list[i], stype='S'))

                    # Cells process
                    bounds = ''
                    for num in surf_num_list:
                        bounds = bounds + '&' + str(num)
                    cells = []
                    # the interval cell
                    cells.append(Cell(number=cell_num_list[14], bounds=bounds[1:len(bounds)],
                                      fill=univ.lattice.fill_interval, volume=None, density=None))
                    # the inner cells
                    move_list = [[0, 0, 0], [l, 0, 0], [l, l, 0], [0, l, 0], [0, 0, l], [l, 0, l], [l, l, l], [0, l, l],
                                 [l / 2, 0, l / 2], [l, l / 2, l / 2], [l / 2, l, l / 2], [0, l / 2, l / 2],
                                 [l / 2, l / 2, l], [l / 2, l / 2, 0]]

                    for i in range(14):
                        cells.append(
                            Cell(number=cell_num_list[i], bounds=str(-surf_num_list[0]), fill=univ.lattice.fill_ball,
                                 transformation=Transformation(move=move_list[i]), volume=None, density=None))

                    post_univ.add_cells(cells)
                    self.universes.append(post_univ)
                    # universe 添加完毕

                    univ.lattice = post_lat
                    # lat=5 解析完毕
        return surfs

    def get_max_uni_id(self):
        # get the maximum user id for univs
        uni_max_id = int(0)
        for univ in self.universes:
            if univ.number > uni_max_id:
                uni_max_id = univ.number
        return uni_max_id

    def get_max_cell_id(self):
        # get the maximum user id for cells
        cell_max_id = int(0)
        for univ in self.universes:
            for cell in univ.cells:
                if cell.number > cell_max_id:
                    cell_max_id = cell.number
        return cell_max_id

    def delete_univ(self, uid):
        for univ in self.universes:
            if univ.number == uid:
                self.universes.remove(univ)

    def get_new_univ_id(self, base_id=0):
        if base_id < 0:
            warnings.warn("The ID of Universe should not be negative, use 0 as the base id.")
            base_id = 0
        while base_id in self.univ_dict:
            base_id += 1
        return base_id

    def get_new_cell_id(self, base_id=1):
        if base_id < 1:
            warnings.warn("The ID of Cell should not be smaller than 1, use 1 as the base id.")
            base_id = 1
        while base_id in self.cell_dict:
            base_id += 1
        return base_id

    def __str__(self):
        s = ''
        for univ in self.universes:
            s += str(univ)
        return s

    def __iter__(self):
        for univ in self.universes:
            yield univ





if __name__ == '__main__':
    px = PX(stype='PX', parameters=[1], boundary=1, rotate=[1, 0, 0, 0, 1, 0, 0, 0, 1])
    p = px.transfer()
