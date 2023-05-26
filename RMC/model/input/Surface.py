# -*- coding:utf-8 -*-
# author: LuoHao, Jiangjiaruo, Gouyuanhao
# date: 2021-12-21
import numpy as np
import re
import sys

from RMC.model.input.base import YMLModelObject as BaseModel
import math
from abc import abstractmethod


class Surfaces(BaseModel):
    """
    Attributes
    -----------
    surfaces: list of Surface object
        all surfaces
    """
    yaml_tag = u'!surfaces'

    def __init__(self, surfaces=None):
        self.surfaces = surfaces
        if self.surfaces is None:
            self.surfaces = []

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

    def get_max_surf_id(self):
        """ get the maximum user id for surfs

        Returns
        -------
        maxmium surface id
        """
        surf_max_id = int(0)
        for surf in self.surfaces:
            if surf.number > surf_max_id:
                surf_max_id = surf.number
        return surf_max_id

    def __str__(self):
        s = 'SURFACE\n'
        for surf in self.surfaces:
            s += str(surf)
        s += '\n\n'
        return s


class Surface(BaseModel):
    """
    Attributes
    -----------
    number: int
        surface id
    type: string
        surface type
    parameters: list of float
        surface parameters
    boundary: int
        surface boundary condition
    pair: int
        pair surface id
    value: list of float
        value card
    time: list of float
        time card
    rotate: list of float
        surface rotate matrix
    rotate_angle: list of float
        surface rotate angle
    move: list of float
        surface translation parameters

    """
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
        'PAIR': [int],
        'VALUE': ['list', float, -1],
        'TIME': ['list', float, -1],
        'ROTATE': ['list', float, 9],
        'ROTATEANGLE': ['list', float, 3],
        'MOVE': ['list', float, 3]
    }

    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        self._number = number
        self._type = stype
        self._parameters = parameters
        self._boundary = boundary
        self._pair = pair
        self._value = value
        self._time = time
        self._rotate = rotate
        self._rotate_angle = rotate_angle
        self._move = move

    @staticmethod
    def externalization(number=None, type=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                        rotate=None, rotate_angle=None, move=None):
        """ externalize surface object according to surface type

        Returns
        -------
        the externalized surface object
        """
        attribute = re.sub('/', '_', type)
        return getattr(sys.modules[__name__], attribute)(number=number, stype=type, parameters=parameters,
                                                         boundary=boundary, pair=pair, time=time,
                                                         value=value, rotate=rotate,
                                                         rotate_angle=rotate_angle, move=move)

    @property
    def number(self):
        return self._number

    @number.setter
    def number(self, number):
        self._number = number

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, stype):
        self._type = stype

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, para):
        self._parameters = para

    @property
    def boundary(self):
        return self._boundary

    @boundary.setter
    def boundary(self, boundary):
        self._boundary = boundary

    @property
    def pair(self):
        return self._pair

    @pair.setter
    def pair(self, pair):
        self._pair = pair

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, v):
        self._value = v

    @property
    def rotate(self):
        return self._rotate

    @rotate.setter
    def rotate(self, r):
        self._rotate = r

    @property
    def rotate_angle(self):
        return self._rotate_angle

    @rotate_angle.setter
    def rotate_angle(self, r):
        self._rotate_angle = r

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, t):
        self._time = t

    @property
    def move(self):
        return self._move

    @move.setter
    def move(self, m):
        self._move = m

    def check(self):
        assert self._number >= 0

        if self._rotate is not None and self._rotate_angle is not None:
            raise ValueError("the rotate and rotate_angle card can not be provided simutanuously!\n")

    def __str__(self):
        card = 'Surf ' + str(self._number) + ' ' + self._type.lower() + ' '
        surf_para = ' '.join([f'{x:.12g}' for x in self._parameters])
        card += surf_para
        if self._boundary is not None:
            card += ' bc = ' + str(self._boundary)
        if self._pair is not None:
            card += ' pair = ' + str(self._pair)
        if self._value is not None:
            card += ' value = ' + ' '.join([str(x) for x in self._value])
            card += ' time = ' + ' '.join([str(x) for x in self._time])
        card += '\n'
        return card

    @staticmethod
    def surface_list(surf_list, macrobodies):
        """ modify the surface id in macrobodies to one in surfaces

        Parameters
        ----------
        surf_list: surface ids involving id using the macrobody format
        macrobodies: dict of macrobodie object

        Returns
        -------
        modified surface lists

        """
        for i, surface in enumerate(surf_list):
            try:
                float(surface)
                if abs(surface) - int(abs(surface)) == 0:
                    if int(abs(surface)) in macrobodies.keys():
                        surf_list[i] = macrobodies[int(abs(surface))].surf_ids[0]
                    else:
                        continue
                else:
                    surf_list[i] = macrobodies[int(abs(surface))].surf_ids[int(str(surface).split('.')[1]) - 1]
                if surface < 0:
                    surf_list[i] = -surf_list[i]
            except ValueError:
                continue
        return surf_list

    @staticmethod
    def rotate_matrix(alpha_degree, beta_degree, gamma_degree):
        """ generate rotation matrix accorting rotation matrix.
            In implementation of RMC, the extrinsic rotation is adopted and
            the rotation order is z-y-x()
            reference: https://en.wikipedia.org/wiki/Euler_angles

        Parameters
        ----------
        alpha_degree: procession angle
        beta_degree: nutation angle
        gamma_degree：spin angle

        Returns
        -------
        rotation matrix
        """
        rotatematrix = np.zeros((3, 3))

        alpha = alpha_degree * math.pi / 180
        beta = beta_degree * math.pi / 180
        gamma = gamma_degree * math.pi / 180

        rotatematrix[0][0] = math.cos(alpha) * math.cos(beta)
        rotatematrix[0][1] = math.cos(alpha) * math.sin(beta) * math.sin(gamma) - math.sin(alpha) * math.cos(gamma)
        rotatematrix[0][2] = math.cos(alpha) * math.sin(beta) * math.cos(gamma) + math.sin(alpha) * math.sin(gamma)
        rotatematrix[1][0] = math.sin(alpha) * math.cos(beta)
        rotatematrix[1][1] = math.sin(alpha) * math.sin(beta) * math.sin(gamma) + math.cos(alpha) * math.cos(gamma)
        rotatematrix[1][2] = math.sin(alpha) * math.sin(beta) * math.cos(gamma) - math.cos(alpha) * math.sin(gamma)
        rotatematrix[2][0] = -math.sin(beta)
        rotatematrix[2][1] = math.cos(beta) * math.sin(gamma)
        rotatematrix[2][2] = math.cos(beta) * math.cos(gamma)
        return rotatematrix

    def transfer(self):
        """
        Returns
        -------
        Surface object after rotation and translation
        """
        surface = self.rotation()
        return surface.translation()

    @abstractmethod
    def translation(self):
        pass

    @abstractmethod
    def rotation(self):
        pass


class P(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        parameters = [0] * 4
        parameters[0] = self.parameters[0] * rotatematrix[0][0] + self.parameters[1] * rotatematrix[1][0] + \
                        self.parameters[2] * rotatematrix[2][0]
        parameters[1] = self.parameters[0] * rotatematrix[0][1] + self.parameters[1] * rotatematrix[1][1] + \
                        self.parameters[2] * rotatematrix[2][1]
        parameters[2] = self.parameters[0] * rotatematrix[0][2] + self.parameters[1] * rotatematrix[1][2] + \
                        self.parameters[2] * rotatematrix[2][2]
        parameters[3] = self.parameters[3]

        return P(number=self.number, stype='P', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[3] = self.parameters[0] * self.move[0] + self.parameters[1] * self.move[1] + \
                                 self.parameters[2] * self.move[2]
            self.move = None
        return self


class PX(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        type = 'P'
        parameters = [0] * 4
        parameters[0] = rotatematrix[0][0]
        parameters[1] = rotatematrix[0][1]
        parameters[2] = rotatematrix[0][2]
        parameters[3] = self.parameters[0]

        return P(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.move = None
        return self


class PY(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'P'

        parameters = [0] * 4
        parameters[0] = rotatematrix[1][0]
        parameters[1] = rotatematrix[1][1]
        parameters[2] = rotatematrix[1][2]
        parameters[3] = self.parameters[0]
        return P(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[1]
            self.move = None
        return self


class PZ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'P'

        parameters = [0] * 4
        parameters[0] = rotatematrix[2][0]
        parameters[1] = rotatematrix[2][1]
        parameters[2] = rotatematrix[2][2]
        parameters[3] = self.parameters[0]
        return P(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[2]
            self.move = None
        return self


class SO(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        return self

    def translation(self):
        if self.move is None:
            return self

        parameters = [0] * 4

        parameters[0] = self.move[0]
        parameters[1] = self.move[1]
        parameters[2] = self.move[2]
        parameters[3] = self.parameters[0]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time)


class S(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        parameters = [0] * 4

        parameters[0] = self.parameters[0] * rotatematrix[0][0] + self.parameters[1] * rotatematrix[1][0] + \
                        self.parameters[2] * rotatematrix[2][0]
        parameters[1] = self.parameters[0] * rotatematrix[0][1] + self.parameters[1] * rotatematrix[1][1] + \
                        self.parameters[2] * rotatematrix[2][1]
        parameters[2] = self.parameters[0] * rotatematrix[0][2] + self.parameters[1] * rotatematrix[1][2] + \
                        self.parameters[2] * rotatematrix[2][2]
        parameters[3] = self.parameters[3]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] = self.parameters[0] + self.move[0]
            self.parameters[1] = self.parameters[1] + self.move[1]
            self.parameters[2] = self.parameters[2] + self.move[2]
            self.move = None
        return self


class SX(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        parameters = [0] * 4
        parameters[0] = self.parameters[0] * rotatematrix[0][0]
        parameters[1] = self.parameters[0] * rotatematrix[0][1]
        parameters[2] = self.parameters[0] * rotatematrix[0][2]
        parameters[3] = self.parameters[1]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self

        parameters = [0] * 4
        parameters[0] = self.parameters[0] + self.move[0]
        parameters[1] = self.move[1]
        parameters[2] = self.move[2]
        parameters[3] = self.parameters[1]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time)


class SY(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        parameters = [0] * 4

        parameters[0] = self.parameters[0] * rotatematrix[1][0]
        parameters[1] = self.parameters[0] * rotatematrix[1][1]
        parameters[2] = self.parameters[0] * rotatematrix[1][2]
        parameters[3] = self.parameters[1]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self

        parameters = [0] * 4
        parameters[0] = self.move[0]
        parameters[1] = self.parameters[0] + self.move[1]
        parameters[2] = self.move[2]
        parameters[3] = self.parameters[1]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time)


class SZ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        parameters = [0] * 4

        parameters[0] = self.parameters[0] * rotatematrix[2][0]
        parameters[1] = self.parameters[0] * rotatematrix[2][1]
        parameters[2] = self.parameters[0] * rotatematrix[2][2]
        parameters[3] = self.parameters[1]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self

        parameters = [0] * 4
        parameters[0] = self.move[0]
        parameters[1] = self.move[1]
        parameters[2] = self.parameters[0] + self.move[2]
        parameters[3] = self.parameters[1]

        return S(number=self.number, stype='S', parameters=parameters, boundary=self.boundary, pair=self.pair,
                 value=self.value, time=self.time)


class C_X(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = 1 - rotatematrix[0][0] * rotatematrix[0][0]
        parameters[1] = 1 - rotatematrix[0][1] * rotatematrix[0][1]
        parameters[2] = 1 - rotatematrix[0][2] * rotatematrix[0][2]
        parameters[3] = (rotatematrix[1][0] * rotatematrix[1][1] + rotatematrix[2][0] * rotatematrix[2][1]) * 2
        parameters[4] = (rotatematrix[1][1] * rotatematrix[1][2] + rotatematrix[2][1] * rotatematrix[2][2]) * 2
        parameters[5] = (rotatematrix[1][0] * rotatematrix[1][2] + rotatematrix[2][0] * rotatematrix[2][2]) * 2
        parameters[6] = -(rotatematrix[1][0] * self.parameters[0] + rotatematrix[2][0] * self.parameters[1]) * 2
        parameters[7] = -(rotatematrix[1][1] * self.parameters[0] + rotatematrix[2][1] * self.parameters[1]) * 2
        parameters[8] = -(rotatematrix[1][2] * self.parameters[0] + rotatematrix[2][2] * self.parameters[1]) * 2
        parameters[9] = -(self.parameters[2] * self.parameters[2] - self.parameters[0] * self.parameters[0] - \
                          self.parameters[1] * self.parameters[1])

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[1]
            self.parameters[1] += self.move[2]
            self.move = None
        return self


class C_Y(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'

        parameters = [0] * 10

        parameters[0] = 1 - rotatematrix[1][0] * rotatematrix[1][0]
        parameters[1] = 1 - rotatematrix[1][1] * rotatematrix[1][1]
        parameters[2] = 1 - rotatematrix[1][2] * rotatematrix[1][2]
        parameters[3] = (rotatematrix[0][0] * rotatematrix[0][1] + rotatematrix[2][0] * rotatematrix[2][1]) * 2
        parameters[4] = (rotatematrix[0][1] * rotatematrix[0][2] + rotatematrix[2][1] * rotatematrix[2][2]) * 2
        parameters[5] = (rotatematrix[0][0] * rotatematrix[0][2] + rotatematrix[2][0] * rotatematrix[2][2]) * 2
        parameters[6] = -(rotatematrix[0][0] * self.parameters[0] + rotatematrix[2][0] * self.parameters[1]) * 2
        parameters[7] = -(rotatematrix[0][1] * self.parameters[0] + rotatematrix[2][1] * self.parameters[1]) * 2
        parameters[8] = -(rotatematrix[0][2] * self.parameters[0] + rotatematrix[2][2] * self.parameters[1]) * 2
        parameters[9] = -(self.parameters[2] * self.parameters[2] - self.parameters[0] * self.parameters[0] - \
                          self.parameters[1] * self.parameters[1])

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[2]
            self.move = None
        return self


class C_Z(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = 1 - rotatematrix[2][0] * rotatematrix[2][0]
        parameters[1] = 1 - rotatematrix[2][1] * rotatematrix[2][1]
        parameters[2] = 1 - rotatematrix[2][2] * rotatematrix[2][2]
        parameters[3] = (rotatematrix[1][0] * rotatematrix[1][1] + rotatematrix[0][0] * rotatematrix[0][1]) * 2
        parameters[4] = (rotatematrix[1][1] * rotatematrix[1][2] + rotatematrix[0][1] * rotatematrix[0][2]) * 2
        parameters[5] = (rotatematrix[1][0] * rotatematrix[1][2] + rotatematrix[0][0] * rotatematrix[0][2]) * 2
        parameters[6] = -(rotatematrix[0][0] * self.parameters[0] + rotatematrix[1][0] * self.parameters[1]) * 2
        parameters[7] = -(rotatematrix[0][1] * self.parameters[0] + rotatematrix[1][1] * self.parameters[1]) * 2
        parameters[8] = -(rotatematrix[0][2] * self.parameters[0] + rotatematrix[1][2] * self.parameters[1]) * 2
        parameters[9] = -(self.parameters[2] * self.parameters[2] - self.parameters[0] * self.parameters[0] - \
                          self.parameters[1] * self.parameters[1])

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.move = None
        return self


class CX(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = 1 - rotatematrix[0][0] * rotatematrix[0][0]
        parameters[1] = 1 - rotatematrix[0][1] * rotatematrix[0][1]
        parameters[2] = 1 - rotatematrix[0][2] * rotatematrix[0][2]
        parameters[3] = (rotatematrix[1][0] * rotatematrix[1][1] + rotatematrix[2][0] * rotatematrix[2][1]) * 2
        parameters[4] = (rotatematrix[1][1] * rotatematrix[1][2] + rotatematrix[2][1] * rotatematrix[2][2]) * 2
        parameters[5] = (rotatematrix[1][0] * rotatematrix[1][2] + rotatematrix[2][0] * rotatematrix[2][2]) * 2
        parameters[6] = 0
        parameters[7] = 0
        parameters[8] = 0
        parameters[9] = -self.parameters[0] ** 2

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self
        parameters = [0] * 3
        parameters[0] = self.move[1]
        parameters[1] = self.move[2]
        parameters[2] = self.parameters[0]
        return C_X(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time)


class CY(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = 1 - rotatematrix[1][0] * rotatematrix[1][0]
        parameters[1] = 1 - rotatematrix[1][1] * rotatematrix[1][1]
        parameters[2] = 1 - rotatematrix[1][2] * rotatematrix[1][2]
        parameters[3] = (rotatematrix[0][0] * rotatematrix[0][1] + rotatematrix[2][0] * rotatematrix[2][1]) * 2
        parameters[4] = (rotatematrix[0][1] * rotatematrix[0][2] + rotatematrix[2][1] * rotatematrix[2][2]) * 2
        parameters[5] = (rotatematrix[0][0] * rotatematrix[0][2] + rotatematrix[2][0] * rotatematrix[2][2]) * 2
        parameters[6] = 0
        parameters[7] = 0
        parameters[8] = 0
        parameters[9] = -self.parameters[0] ** 2

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self
        parameters = [0] * 3
        parameters[0] = self.move[0]
        parameters[1] = self.move[2]
        parameters[2] = self.parameters[0]
        return C_Y(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time)


class CZ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = 1 - rotatematrix[2][0] * rotatematrix[2][0]
        parameters[1] = 1 - rotatematrix[2][1] * rotatematrix[2][1]
        parameters[2] = 1 - rotatematrix[2][2] * rotatematrix[2][2]
        parameters[3] = (rotatematrix[1][0] * rotatematrix[1][1] + rotatematrix[0][0] * rotatematrix[0][1]) * 2
        parameters[4] = (rotatematrix[1][1] * rotatematrix[1][2] + rotatematrix[0][1] * rotatematrix[0][2]) * 2
        parameters[5] = (rotatematrix[1][0] * rotatematrix[1][2] + rotatematrix[0][0] * rotatematrix[0][2]) * 2
        parameters[6] = 0
        parameters[7] = 0
        parameters[8] = 0
        parameters[9] = -self.parameters[0] ** 2

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self
        parameters = [0] * 3
        parameters[0] = self.move[0]
        parameters[1] = self.move[1]
        parameters[2] = self.parameters[0]
        return C_Z(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time)


class K_X(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'KGQ'
        parameters = [0] * 16

        parameters[0] = 1 - (self.parameters[3] + 1) * rotatematrix[0][0] * rotatematrix[0][0]
        parameters[1] = 1 - (self.parameters[3] + 1) * rotatematrix[0][1] * rotatematrix[0][1]
        parameters[2] = 1 - (self.parameters[3] + 1) * rotatematrix[0][2] * rotatematrix[0][2]
        parameters[3] = -2 * (
                self.parameters[3] * rotatematrix[0][0] * rotatematrix[0][1] - rotatematrix[1][0] * rotatematrix[1][
            1] - rotatematrix[2][0] * rotatematrix[2][1])
        parameters[4] = -2 * (
                self.parameters[3] * rotatematrix[0][1] * rotatematrix[0][2] - rotatematrix[1][1] * rotatematrix[1][
            2] - rotatematrix[2][1] * rotatematrix[2][2])
        parameters[5] = -2 * (
                self.parameters[3] * rotatematrix[0][0] * rotatematrix[0][2] - rotatematrix[1][0] * rotatematrix[1][
            2] - rotatematrix[2][0] * rotatematrix[2][2])
        parameters[6] = -2 * (-self.parameters[3] * rotatematrix[0][0] * self.parameters[0] + rotatematrix[1][0] *
                              self.parameters[1] + rotatematrix[2][0] * self.parameters[2])
        parameters[7] = -2 * (-self.parameters[3] * rotatematrix[0][1] * self.parameters[0] + rotatematrix[1][1] *
                              self.parameters[1] + rotatematrix[2][1] * self.parameters[2])
        parameters[8] = -2 * (-self.parameters[3] * rotatematrix[0][2] * self.parameters[0] + rotatematrix[1][2] *
                              self.parameters[1] + rotatematrix[2][2] * self.parameters[2])
        parameters[9] = -(
                self.parameters[3] * self.parameters[0] ** 2 - self.parameters[1] ** 2 - self.parameters[2] ** 2)

        parameters[10] = rotatematrix[0][0] * self.parameters[0] + rotatematrix[1][0] * self.parameters[1] + \
                         rotatematrix[2][0] * self.parameters[2]
        parameters[11] = rotatematrix[0][1] * self.parameters[0] + rotatematrix[1][1] * self.parameters[1] + \
                         rotatematrix[2][1] * self.parameters[2]
        parameters[12] = rotatematrix[0][2] * self.parameters[0] + rotatematrix[1][2] * self.parameters[1] + \
                         rotatematrix[2][2] * self.parameters[2]
        parameters[13] = rotatematrix[0][0] * self.parameters[4]
        parameters[14] = rotatematrix[0][1] * self.parameters[4]
        parameters[15] = rotatematrix[0][2] * self.parameters[4]

        return KGQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.parameters[2] += self.move[2]
            self.move = None
        return self


class K_Y(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        type = 'KGQ'

        parameters = [0] * 16

        parameters[0] = 1 - (self.parameters[3] + 1) * rotatematrix[1][0] * rotatematrix[1][0]
        parameters[1] = 1 - (self.parameters[3] + 1) * rotatematrix[1][1] * rotatematrix[1][1]
        parameters[2] = 1 - (self.parameters[3] + 1) * rotatematrix[1][2] * rotatematrix[1][2]
        parameters[3] = -2 * (
                self.parameters[3] * rotatematrix[1][0] * rotatematrix[1][1] - rotatematrix[0][0] * rotatematrix[0][
            1] - rotatematrix[2][0] * rotatematrix[2][1])
        parameters[4] = -2 * (
                self.parameters[3] * rotatematrix[1][1] * rotatematrix[1][2] - rotatematrix[0][1] * rotatematrix[0][
            2] - rotatematrix[2][1] * rotatematrix[2][2])
        parameters[5] = -2 * (
                self.parameters[3] * rotatematrix[1][0] * rotatematrix[1][2] - rotatematrix[0][0] * rotatematrix[0][
            2] - rotatematrix[2][0] * rotatematrix[2][2])
        parameters[6] = -2 * (-self.parameters[3] * rotatematrix[1][0] * self.parameters[1] + rotatematrix[0][0] *
                              self.parameters[0] + rotatematrix[2][0] * self.parameters[2])
        parameters[7] = -2 * (-self.parameters[3] * rotatematrix[1][1] * self.parameters[1] + rotatematrix[0][1] *
                              self.parameters[0] + rotatematrix[2][1] * self.parameters[2])
        parameters[8] = -2 * (-self.parameters[3] * rotatematrix[1][2] * self.parameters[1] + rotatematrix[0][2] *
                              self.parameters[0] + rotatematrix[2][2] * self.parameters[2])
        parameters[9] = -(
                self.parameters[3] * self.parameters[1] ** 2 - self.parameters[0] ** 2 - self.parameters[2] ** 2)
        parameters[10] = rotatematrix[0][0] * self.parameters[0] + rotatematrix[1][0] * self.parameters[1] + \
                         rotatematrix[2][0] * self.parameters[2]
        parameters[11] = rotatematrix[0][1] * self.parameters[0] + rotatematrix[1][1] * self.parameters[1] + \
                         rotatematrix[2][1] * self.parameters[2]
        parameters[12] = rotatematrix[0][2] * self.parameters[0] + rotatematrix[1][2] * self.parameters[1] + \
                         rotatematrix[2][2] * self.parameters[2]
        parameters[13] = rotatematrix[1][0] * self.parameters[4]
        parameters[14] = rotatematrix[1][1] * self.parameters[4]
        parameters[15] = rotatematrix[1][2] * self.parameters[4]

        return KGQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.parameters[2] += self.move[2]
            self.move = None
        return self


class K_Z(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        type = 'KGQ'

        parameters = [0] * 16

        parameters[0] = 1 - (self.parameters[3] + 1) * rotatematrix[2][0] * rotatematrix[2][0]
        parameters[1] = 1 - (self.parameters[3] + 1) * rotatematrix[2][1] * rotatematrix[2][1]
        parameters[2] = 1 - (self.parameters[3] + 1) * rotatematrix[2][2] * rotatematrix[2][2]
        parameters[3] = -2 * (
                self.parameters[3] * rotatematrix[2][0] * rotatematrix[2][1] - rotatematrix[0][0] * rotatematrix[0][
            1] - rotatematrix[1][0] * rotatematrix[1][1])
        parameters[4] = -2 * (
                self.parameters[3] * rotatematrix[2][1] * rotatematrix[2][2] - rotatematrix[0][1] * rotatematrix[0][
            2] - rotatematrix[1][1] * rotatematrix[1][2])
        parameters[5] = -2 * (
                self.parameters[3] * rotatematrix[2][0] * rotatematrix[2][2] - rotatematrix[0][0] * rotatematrix[0][
            2] - rotatematrix[1][0] * rotatematrix[1][2])
        parameters[6] = -2 * (-self.parameters[3] * rotatematrix[2][0] * self.parameters[2] + rotatematrix[0][0] *
                              self.parameters[0] + rotatematrix[1][0] * self.parameters[1])
        parameters[7] = -2 * (-self.parameters[3] * rotatematrix[2][1] * self.parameters[2] + rotatematrix[0][1] *
                              self.parameters[0] + rotatematrix[1][1] * self.parameters[1])
        parameters[8] = -2 * (-self.parameters[3] * rotatematrix[2][2] * self.parameters[2] + rotatematrix[0][2] *
                              self.parameters[0] + rotatematrix[1][2] * self.parameters[1])
        parameters[9] = -(
                self.parameters[3] * self.parameters[2] ** 2 - self.parameters[0] ** 2 - self.parameters[1] ** 2)
        parameters[10] = rotatematrix[0][0] * self.parameters[0] + rotatematrix[1][0] * self.parameters[1] + \
                         rotatematrix[2][0] * self.parameters[2]
        parameters[11] = rotatematrix[0][1] * self.parameters[0] + rotatematrix[1][1] * self.parameters[1] + \
                         rotatematrix[2][1] * self.parameters[2]
        parameters[12] = rotatematrix[0][2] * self.parameters[0] + rotatematrix[1][2] * self.parameters[1] + \
                         rotatematrix[2][2] * self.parameters[2]
        parameters[13] = rotatematrix[2][0] * self.parameters[4]
        parameters[14] = rotatematrix[2][1] * self.parameters[4]
        parameters[15] = rotatematrix[2][2] * self.parameters[4]

        return KGQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.parameters[2] += self.move[2]
            self.move = None
        return self


class KX(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])
        type = 'KGQ'
        parameters = [0] * 16

        parameters[0] = 1 - (self.parameters[1] + 1) * rotatematrix[0][0] * rotatematrix[0][0]
        parameters[1] = 1 - (self.parameters[1] + 1) * rotatematrix[0][1] * rotatematrix[0][1]
        parameters[2] = 1 - (self.parameters[1] + 1) * rotatematrix[0][2] * rotatematrix[0][2]
        parameters[3] = -2 * (
                self.parameters[1] * rotatematrix[0][0] * rotatematrix[0][1] - rotatematrix[1][0] * rotatematrix[1][
            1] - rotatematrix[2][0] * rotatematrix[2][1])
        parameters[4] = -2 * (
                self.parameters[1] * rotatematrix[0][1] * rotatematrix[0][2] - rotatematrix[1][1] * rotatematrix[1][
            2] - rotatematrix[2][1] * rotatematrix[2][2])
        parameters[5] = -2 * (
                self.parameters[1] * rotatematrix[0][0] * rotatematrix[0][2] - rotatematrix[1][0] * rotatematrix[1][
            2] - rotatematrix[2][0] * rotatematrix[2][2])
        parameters[6] = -2 * (-self.parameters[1] * rotatematrix[0][0] * self.parameters[0])
        parameters[7] = -2 * (-self.parameters[1] * rotatematrix[0][1] * self.parameters[0])
        parameters[8] = -2 * (-self.parameters[1] * rotatematrix[0][2] * self.parameters[0])
        parameters[9] = -self.parameters[1] * self.parameters[0] ** 2

        parameters[10] = rotatematrix[0][0] * self.parameters[0]
        parameters[11] = rotatematrix[0][1] * self.parameters[0]
        parameters[12] = rotatematrix[0][2] * self.parameters[0]
        parameters[13] = rotatematrix[0][0] * self.parameters[2]
        parameters[14] = rotatematrix[0][1] * self.parameters[2]
        parameters[15] = rotatematrix[0][2] * self.parameters[2]

        return KGQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self
        parameters = [0] * 5

        parameters[0] = self.parameters[0] + self.move[0]
        parameters[1] = self.move[1]
        parameters[2] = self.move[2]
        parameters[3] = self.parameters[1]
        parameters[4] = self.parameters[2]

        return K_X(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time)


class KY(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'KGQ'

        parameters = [0] * 16

        parameters[0] = 1 - (self.parameters[1] + 1) * rotatematrix[1][0] * rotatematrix[1][0]
        parameters[1] = 1 - (self.parameters[1] + 1) * rotatematrix[1][1] * rotatematrix[1][1]
        parameters[2] = 1 - (self.parameters[1] + 1) * rotatematrix[1][2] * rotatematrix[1][2]
        parameters[3] = -2 * (
                self.parameters[1] * rotatematrix[1][0] * rotatematrix[1][1] - rotatematrix[0][0] * rotatematrix[0][
            1] - rotatematrix[2][0] * rotatematrix[2][1])
        parameters[4] = -2 * (
                self.parameters[1] * rotatematrix[1][1] * rotatematrix[1][2] - rotatematrix[0][1] * rotatematrix[0][
            2] - rotatematrix[2][1] * rotatematrix[2][2])
        parameters[5] = -2 * (
                self.parameters[1] * rotatematrix[1][0] * rotatematrix[1][2] - rotatematrix[0][0] * rotatematrix[0][
            2] - rotatematrix[2][0] * rotatematrix[2][2])
        parameters[6] = -2 * (-self.parameters[1] * rotatematrix[1][0] * self.parameters[0])
        parameters[7] = -2 * (-self.parameters[1] * rotatematrix[1][1] * self.parameters[0])
        parameters[8] = -2 * (-self.parameters[1] * rotatematrix[1][2] * self.parameters[0])
        parameters[9] = -(self.parameters[1] * self.parameters[0] ** 2)

        parameters[10] = rotatematrix[1][0] * self.parameters[0]
        parameters[11] = rotatematrix[1][1] * self.parameters[0]
        parameters[12] = rotatematrix[1][2] * self.parameters[0]
        parameters[13] = rotatematrix[1][0] * self.parameters[2]
        parameters[14] = rotatematrix[1][1] * self.parameters[2]
        parameters[15] = rotatematrix[1][2] * self.parameters[2]

        return KGQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self
        parameters = [0] * 5

        parameters[0] = self.move[0]
        parameters[1] = self.parameters[0] + self.move[1]
        parameters[2] = self.move[2]
        parameters[3] = self.parameters[1]
        parameters[4] = self.parameters[2]

        return K_Y(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time)


class KZ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'KGQ'

        parameters = [0] * 16

        parameters[0] = 1 - (self.parameters[1] + 1) * rotatematrix[2][0] * rotatematrix[2][0]
        parameters[1] = 1 - (self.parameters[1] + 1) * rotatematrix[2][1] * rotatematrix[2][1]
        parameters[2] = 1 - (self.parameters[1] + 1) * rotatematrix[2][2] * rotatematrix[2][2]
        parameters[3] = -2 * (
                self.parameters[1] * rotatematrix[2][0] * rotatematrix[2][1] - rotatematrix[0][0] * rotatematrix[0][
            1] - rotatematrix[1][0] * rotatematrix[1][1])
        parameters[4] = -2 * (
                self.parameters[1] * rotatematrix[2][1] * rotatematrix[2][2] - rotatematrix[0][1] * rotatematrix[0][
            2] - rotatematrix[1][1] * rotatematrix[1][2])
        parameters[5] = -2 * (
                self.parameters[1] * rotatematrix[2][0] * rotatematrix[2][2] - rotatematrix[0][0] * rotatematrix[0][
            2] - rotatematrix[1][0] * rotatematrix[1][2])
        parameters[6] = -2 * (-self.parameters[1] * rotatematrix[2][0] * self.parameters[0])
        parameters[7] = -2 * (-self.parameters[1] * rotatematrix[2][1] * self.parameters[0])
        parameters[8] = -2 * (-self.parameters[1] * rotatematrix[2][2] * self.parameters[0])
        parameters[9] = -(self.parameters[1] * self.parameters[0] ** 2)

        parameters[10] = rotatematrix[2][0] * self.parameters[0]
        parameters[11] = rotatematrix[2][1] * self.parameters[0]
        parameters[12] = rotatematrix[2][2] * self.parameters[0]
        parameters[13] = rotatematrix[2][0] * self.parameters[2]
        parameters[14] = rotatematrix[2][1] * self.parameters[2]
        parameters[15] = rotatematrix[2][2] * self.parameters[2]

        return KGQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is None:
            return self
        parameters = [0] * 5

        parameters[0] = self.move[0]
        parameters[1] = self.move[1]
        parameters[2] = self.parameters[0] + self.move[2]
        parameters[3] = self.parameters[1]
        parameters[4] = self.parameters[2]

        return K_Z(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                   value=self.value, time=self.time)


class KGQ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        return self

    def translation(self):
        if self.move is not None:
            self.parameters[9] += self.parameters[0] * self.move[0] * self.move[0] + self.parameters[1] * self.move[1] * \
                                  self.move[1] + self.parameters[2] * self.move[2] * self.move[2] + self.parameters[2] * \
                                  self.move[0] * self.move[1] + self.parameters[4] * self.move[1] * self.move[2] + \
                                  self.parameters[5] * self.move[0] * self.move[2] - self.parameters[6] * self.move[0] - \
                                  self.parameters[7] * self.move[1] - self.parameters[8] * self.move[2]
            self.parameters[6] -= 2 * self.parameters[0] * self.move[0] - self.parameters[3] * self.move[1] - \
                                  self.parameters[5] * self.move[2]
            self.parameters[7] -= 2 * self.parameters[1] * self.move[1] - self.parameters[3] * self.move[0] - \
                                  self.parameters[4] * self.move[2]
            self.parameters[8] -= 2 * self.parameters[2] * self.move[2] - self.parameters[4] * self.move[1] - \
                                  self.parameters[5] * self.move[0]
            self.parameters[10] += self.move[0]
            self.parameters[11] += self.move[1]
            self.parameters[12] += self.move[2]
            # 重置move參數
            self.move = None
        return self


class SQ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = self.parameters[0] * rotatematrix[0][0] * rotatematrix[0][0] + self.parameters[1] * \
                        rotatematrix[1][0] * rotatematrix[1][0] + self.parameters[2] * rotatematrix[2][0] * \
                        rotatematrix[2][0]
        parameters[1] = self.parameters[0] * rotatematrix[0][1] * rotatematrix[0][1] + self.parameters[1] * \
                        rotatematrix[1][1] * rotatematrix[1][1] + self.parameters[2] * rotatematrix[2][1] * \
                        rotatematrix[2][1]
        parameters[2] = self.parameters[0] * rotatematrix[0][2] * rotatematrix[0][2] + self.parameters[1] * \
                        rotatematrix[1][2] * rotatematrix[1][2] + self.parameters[2] * rotatematrix[2][2] * \
                        rotatematrix[2][2]
        parameters[3] = 2 * (
                self.parameters[0] * rotatematrix[0][0] * rotatematrix[0][1] + self.parameters[1] * rotatematrix[1][
            0] * rotatematrix[1][1] + self.parameters[2] * rotatematrix[2][0] * rotatematrix[2][1])
        parameters[4] = 2 * (
                self.parameters[0] * rotatematrix[0][1] * rotatematrix[0][2] + self.parameters[1] * rotatematrix[1][
            1] * rotatematrix[1][2] + self.parameters[2] * rotatematrix[2][1] * rotatematrix[2][2])
        parameters[5] = 2 * (
                self.parameters[0] * rotatematrix[0][0] * rotatematrix[0][2] + self.parameters[1] * rotatematrix[1][
            0] * rotatematrix[1][2] + self.parameters[2] * rotatematrix[2][0] * rotatematrix[2][2])
        parameters[6] = 2 * (
                self.parameters[3] * rotatematrix[0][0] + self.parameters[4] * rotatematrix[1][0] + self.parameters[
            5] * rotatematrix[2][0] - self.parameters[0] * rotatematrix[0][0] * self.parameters[7] -
                self.parameters[1] * rotatematrix[1][0] * self.parameters[8] - self.parameters[2] * rotatematrix[2][
                    0] * self.parameters[9])
        parameters[7] = 2 * (
                self.parameters[3] * rotatematrix[0][1] + self.parameters[4] * rotatematrix[1][1] + self.parameters[
            5] * rotatematrix[2][1] - self.parameters[0] * rotatematrix[0][1] * self.parameters[7] -
                self.parameters[1] * rotatematrix[1][1] * self.parameters[8] - self.parameters[2] * rotatematrix[2][
                    1] * self.parameters[9])
        parameters[8] = 2 * (
                self.parameters[3] * rotatematrix[0][2] + self.parameters[4] * rotatematrix[1][2] + self.parameters[
            5] * rotatematrix[2][2] - self.parameters[0] * rotatematrix[0][2] * self.parameters[7] -
                self.parameters[1] * rotatematrix[1][2] * self.parameters[8] - self.parameters[2] * rotatematrix[2][
                    2] * self.parameters[9])
        parameters[9] = self.parameters[0] * self.parameters[7] ** 2 + self.parameters[1] * self.parameters[8] ** 2 + \
                        self.parameters[2] * self.parameters[9] ** 2 - 2 * (
                                self.parameters[3] * self.parameters[7] + self.parameters[4] * self.parameters[8] +
                                self.parameters[5] * self.parameters[9]) + self.parameters[6]

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[7] = self.parameters[7] + self.move[0]
            self.parameters[8] = self.parameters[8] + self.move[1]
            self.parameters[9] = self.parameters[9] + self.move[2]
        return self


class GQ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is None and self.rotate_angle is None:
            return self
        rotatematrix = self.rotate if self.rotate is not None else Surface.rotate_matrix(self.rotate_angle[0],
                                                                                         self.rotate_angle[1],
                                                                                         self.rotate_angle[2])

        type = 'GQ'
        parameters = [0] * 10

        parameters[0] = self.parameters[0] * rotatematrix[0][0] * rotatematrix[0][0] + self.parameters[1] * \
                        rotatematrix[1][0] * rotatematrix[1][0] + self.parameters[2] * rotatematrix[2][0] * \
                        rotatematrix[2][0] + self.parameters[3] * rotatematrix[0][0] * rotatematrix[1][0] + \
                        self.parameters[4] * rotatematrix[1][0] * rotatematrix[2][0] + self.parameters[5] * \
                        rotatematrix[0][0] * rotatematrix[2][0]
        parameters[1] = self.parameters[0] * rotatematrix[0][1] * rotatematrix[0][1] + self.parameters[1] * \
                        rotatematrix[1][1] * rotatematrix[1][1] + self.parameters[2] * rotatematrix[2][1] * \
                        rotatematrix[2][1] + self.parameters[3] * rotatematrix[0][1] * rotatematrix[1][1] + \
                        self.parameters[4] * rotatematrix[1][1] * rotatematrix[2][1] + self.parameters[5] * \
                        rotatematrix[0][1] * rotatematrix[2][1]
        parameters[2] = self.parameters[0] * rotatematrix[0][2] * rotatematrix[0][2] + self.parameters[1] * \
                        rotatematrix[1][2] * rotatematrix[1][2] + self.parameters[2] * rotatematrix[2][2] * \
                        rotatematrix[2][2] + self.parameters[3] * rotatematrix[0][2] * rotatematrix[1][2] + \
                        self.parameters[4] * rotatematrix[1][2] * rotatematrix[2][2] + self.parameters[5] * \
                        rotatematrix[0][2] * rotatematrix[2][2]
        parameters[3] = 2 * (
                self.parameters[0] * rotatematrix[0][0] * rotatematrix[0][1] + self.parameters[1] * rotatematrix[1][
            0] * rotatematrix[1][1] + self.parameters[2] * rotatematrix[2][0] * rotatematrix[2][1]) + \
                        self.parameters[3] * (
                                rotatematrix[0][0] * rotatematrix[1][1] + rotatematrix[0][1] * rotatematrix[1][0]) + \
                        self.parameters[4] * (
                                rotatematrix[1][0] * rotatematrix[2][1] + rotatematrix[1][1] * rotatematrix[2][0]) + \
                        self.parameters[5] * (
                                rotatematrix[0][0] * rotatematrix[2][1] + rotatematrix[0][1] * rotatematrix[2][0])
        parameters[4] = 2 * (
                self.parameters[0] * rotatematrix[0][1] * rotatematrix[0][2] + self.parameters[1] * rotatematrix[1][
            1] * rotatematrix[1][2] + self.parameters[2] * rotatematrix[2][1] * rotatematrix[2][2]) + \
                        self.parameters[3] * (
                                rotatematrix[0][1] * rotatematrix[1][2] + rotatematrix[0][2] * rotatematrix[1][1]) + \
                        self.parameters[4] * (
                                rotatematrix[1][1] * rotatematrix[2][2] + rotatematrix[1][2] * rotatematrix[2][1]) + \
                        self.parameters[5] * (
                                rotatematrix[0][1] * rotatematrix[2][2] + rotatematrix[0][2] * rotatematrix[2][1])
        parameters[5] = 2 * (
                self.parameters[0] * rotatematrix[0][0] * rotatematrix[0][2] + self.parameters[1] * rotatematrix[1][
            0] * rotatematrix[1][2] + self.parameters[2] * rotatematrix[2][0] * rotatematrix[2][2]) + \
                        self.parameters[3] * (
                                rotatematrix[0][0] * rotatematrix[1][2] + rotatematrix[0][2] * rotatematrix[1][0]) + \
                        self.parameters[4] * (
                                rotatematrix[1][0] * rotatematrix[2][2] + rotatematrix[1][2] * rotatematrix[2][0]) + \
                        self.parameters[5] * (
                                rotatematrix[0][0] * rotatematrix[2][2] + rotatematrix[0][2] * rotatematrix[2][0])
        parameters[6] = self.parameters[6] * rotatematrix[0][0] + self.parameters[7] * rotatematrix[1][0] + \
                        self.parameters[8] * rotatematrix[2][0]
        parameters[7] = self.parameters[6] * rotatematrix[0][1] + self.parameters[7] * rotatematrix[1][1] + \
                        self.parameters[8] * rotatematrix[2][1]
        parameters[8] = self.parameters[6] * rotatematrix[0][2] + self.parameters[7] * rotatematrix[1][2] + \
                        self.parameters[8] * rotatematrix[2][2]
        parameters[9] = self.parameters[9]

        return GQ(number=self.number, stype=type, parameters=parameters, boundary=self.boundary, pair=self.pair,
                  value=self.value, time=self.time, move=self.move)

    def translation(self):
        if self.move is not None:
            self.parameters[9] = self.parameters[9] + self.parameters[0] * self.move[0] * self.move[0] + \
                                 self.parameters[1] * self.move[1] * self.move[1] + self.parameters[2] * self.move[2] * \
                                 self.move[2] + self.parameters[3] * self.move[0] * self.move[1] + self.parameters[4] * \
                                 self.move[1] * self.move[2] + self.parameters[5] * self.move[0] * self.move[2] - \
                                 self.parameters[6] * self.move[0] - self.parameters[7] * self.move[1] - \
                                 self.parameters[8] * self.move[2]
            self.parameters[6] = self.parameters[6] - 2 * self.parameters[0] * self.move[0] - self.parameters[3] * \
                                 self.move[1] - \
                                 self.parameters[5] * self.move[2]
            self.parameters[7] = self.parameters[7] - 2 * self.parameters[1] * self.move[1] - self.parameters[3] * \
                                 self.move[0] - \
                                 self.parameters[4] * self.move[2]
            self.parameters[8] = self.parameters[8] - 2 * self.parameters[2] * self.move[2] - self.parameters[4] * \
                                 self.move[1] - \
                                 self.parameters[5] * self.move[0]
            self.move = None
        return self


class TX(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is not None or self.rotate_angle is not None:
            raise NotImplementedError("The rotation of TX surface is not supported!\n")
        return self

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.parameters[2] += self.move[2]
        return self


class TY(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is not None or self.rotate_angle is not None:
            raise NotImplementedError("The rotation of TY surface is not supported!\n")
        return self

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.parameters[2] += self.move[2]
        return self


class TZ(Surface):
    def __init__(self, number=None, stype=None, parameters=None, boundary=None, pair=None, value=None, time=None,
                 rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, stype=stype, parameters=parameters, boundary=boundary, pair=pair, value=value,
                         time=time, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def rotation(self):
        if self.rotate is not None or self.rotate_angle is not None:
            raise NotImplementedError("The rotation of TZ surface is not supported!\n")
        return self

    def translation(self):
        if self.move is not None:
            self.parameters[0] += self.move[0]
            self.parameters[1] += self.move[1]
            self.parameters[2] += self.move[2]
        return self
