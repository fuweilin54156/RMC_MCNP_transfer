from RMC.model.input.Surface import *
import numpy as np
import math
import abc
from RMC.model.input.base import YMLModelObject as BaseModel


class MacroBody(BaseModel):
    yaml_tag = u'!macrobody'

    body_type_para = {
        'RCC': ['list', float, 7],
        'RPP': ['list', float, 6],
        'BOX': ['list', float, 12],
        'SPH': ['list', float, 1 or 4],
        'HEX': ['list', float,  15 or 9],
        'RHP': ['list', float,  9 or 15],
        'REC': ['list', float,  10],
        'TRC': ['list', float, 9],
        'ELL': ['list', float, 7],
        'WED': ['list', float, 12],
        'SEC': ['list', float, 10],
        'TORUS': ['list', float, 8],
        'ROTATE': ['list', float, 9],
        'ROTATEANGLE': ['list', float, 3],
        'MOVE': ['list', float, 3]
    }

    def __init__(self, number=None, type=None, params=None, rotate=None, rotate_angle=None, move=None):
        self.type = type
        self.params = params
        self.body_number = number
        self.surf_ids = []

        self._rotate = rotate
        self._rotate_angle = rotate_angle
        self._move = move

    @property
    def rotate(self):
        return self._rotate

    @property
    def rotate_angle(self):
        return self._rotate_angle

    @property
    def move(self):
        return self._move

    def check(self):
        pass

    @abc.abstractmethod
    # 这个函数用来把各个宏体分解为面，返回面的列表
    def transfer(self, max_surf_id=0):
        pass

    @abc.abstractmethod
    # 返回宏体布尔运算
    def bound(self, in_=True):
        pass

    def __str__(self):
        card = 'Body ' + str(self.body_number) + ' ' + self.type.lower() + ' '
        surf_para = ' '.join([f'{x:.12g}' for x in self.params])
        card += surf_para
        if self._move is not None:
            card += ' move = ' + ' '.join([str(x) for x in self._move])
        if self._rotate is not None:
            card += ' rotate = ' + ' '.join([str(x) for x in self._rotate.reshape([9])])
        card += '\n'
        return card

    @staticmethod
    def externalization(number=None, type=None, parameters=None, rotate=None, rotate_angle=None, move=None):
        """ externalize macrobody object according to macrobody type

        Returns
        -------
        the externalized macrobody object
        """
        attribute = re.sub('/', '_', type, re.I)
        if(attribute == 'RHP'):
            attribute = 'HEX'
        ans = getattr(sys.modules[__name__], attribute, 'None')(number=number, type=type, params=parameters,
                                                         rotate=rotate, rotate_angle=rotate_angle, move=move)

        return ans


class RCC(MacroBody):
    # 圆柱类

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface_model = []
        axial_vector = [self.params[3], self.params[4], self.params[5]]  # 圆柱轴向向量
        co_bottom = [self.params[0], self.params[1], self.params[2]]  # 圆柱底面坐标
        axial_length = np.linalg.norm(axial_vector)  # 圆柱轴向长度
        radius = self.params[6]  # 圆柱半径
        surface = [[] for _ in range(3)]  # 三个面分别的系数列表
        axial_vector = [self.params[3] / axial_length, self.params[4] / axial_length, self.params[5] / axial_length]
        if (axial_vector[0] == 0 or axial_vector[1] == 0) and (axial_vector[0] == 0 or axial_vector[2] == 0) and (
                axial_vector[1] == 0 or axial_vector[2] == 0):
            if [axial_vector[0], axial_vector[1]] == [0, 0]:
                surf_type = ['C/Z', 'P', 'P']
                surface[0].append(co_bottom[0])
                surface[0].append(co_bottom[1])
                surface[0].append(radius)
            elif [axial_vector[0], axial_vector[2]] == [0, 0]:
                surf_type = ['C/Y', 'P', 'P']
                surface[0].append(co_bottom[0])
                surface[0].append(co_bottom[2])
                surface[0].append(radius)
            elif [axial_vector[1], axial_vector[2]] == [0, 0]:
                surf_type = ['C/X', 'P', 'P']
                surface[0].append(co_bottom[1])
                surface[0].append(co_bottom[2])
                surface[0].append(radius)
        else:
            surf_type = ['GQ', 'P', 'P']
            # return gq
            surface[0].append(1 - axial_vector[0] * axial_vector[0])
            surface[0].append(1 - axial_vector[1] * axial_vector[1])
            surface[0].append(1 - axial_vector[2] * axial_vector[2])
            surface[0].append(-2 * axial_vector[0] * axial_vector[1])
            surface[0].append(-2 * axial_vector[1] * axial_vector[2])
            surface[0].append(-2 * axial_vector[0] * axial_vector[2])
            surface[0].append(
                -co_bottom[1] * surface[0][3] - co_bottom[2] * surface[0][5] - 2 * co_bottom[0] * surface[0][0])
            surface[0].append(
                -co_bottom[0] * surface[0][3] - co_bottom[2] * surface[0][4] - 2 * co_bottom[1] * surface[0][1])
            surface[0].append(
                -co_bottom[0] * surface[0][5] - co_bottom[1] * surface[0][4] - 2 * co_bottom[2] * surface[0][2])
            surface[0].append(
                co_bottom[0] * co_bottom[1] * surface[0][3] + co_bottom[1] * co_bottom[2] * surface[0][4] +
                co_bottom[0] * co_bottom[2] * surface[0][5] + co_bottom[0] * co_bottom[0] * surface[0][0] +
                co_bottom[1] * co_bottom[1] * surface[0][1] + co_bottom[2] * co_bottom[2] * surface[0][2] -
                radius * radius)
        # return p
        surface[1].append(axial_vector[0])
        surface[1].append(axial_vector[1])
        surface[1].append(axial_vector[2])
        surface[1].append(axial_vector[0] * (co_bottom[0] + axial_length * axial_vector[0]) + axial_vector[1] * (
                co_bottom[1] + axial_length * axial_vector[1]) + axial_vector[2] * (
                                  co_bottom[2] + axial_length * axial_vector[2]))
        # return p
        surface[2].append(axial_vector[0])
        surface[2].append(axial_vector[1])
        surface[2].append(axial_vector[2])
        surface[2].append(
            axial_vector[0] * co_bottom[0] + axial_vector[1] * co_bottom[1] + axial_vector[2] * co_bottom[2])
        for i in range(0, 3):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{-self.surf_ids[1]}&{self.surf_ids[2]})'
        else:
            return f'({self.surf_ids[0]}:{self.surf_ids[1]}:{-self.surf_ids[2]})'


class RPP(MacroBody):

    # 平行于坐标轴的长方体
    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface = []
        surface_model = []
        surf_type = ['PX', 'PX', 'PY', 'PY', 'PZ', 'PZ']
        for i in range(len(self.params)):
            if i % 2 == 0:
                surface.append(self.params[i + 1])
            else:
                surface.append(self.params[i - 1])
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=[surface[i]],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{self.surf_ids[1]}&{-self.surf_ids[2]}&{self.surf_ids[3]}&{-self.surf_ids[4]}' \
                   f'&{self.surf_ids[5]})'
        else:
            return f'({self.surf_ids[0]}:{-self.surf_ids[1]}:{self.surf_ids[2]}:{-self.surf_ids[3]}:{self.surf_ids[4]}' \
                   f':{-self.surf_ids[5]})'


class BOX(MacroBody):
    # 任意方向的长方体

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface_model = []
        surface = [[] for _ in range(6)]
        box_coord = self.params[0:3]  # 获取其中的一个顶点的坐标
        box_first = np.zeros(3)  # 获取其中的第二个个顶点的坐标
        box_second = np.zeros(3)  # 获取其中的第三个顶点的坐标
        box_third = np.zeros(3)  # 获取其中的第四个个顶点的坐标
        first_vector = self.params[3:6]  # 获取顶点到第一个面的向量
        second_vector = self.params[6:9]  # 获取顶点到第二个面的向量
        third_vector = self.params[9:12]  # 获取顶点到第三个面的向量
        surf_type = ['P', 'P', 'P', 'P', 'P', 'P']
        for i in range(0, 3):
            surface[0].append(first_vector[i])
            surface[1].append(first_vector[i])
            surface[2].append(second_vector[i])
            surface[3].append(second_vector[i])
            surface[4].append(third_vector[i])
            surface[5].append(third_vector[i])
            box_first[i] = box_coord[i] + first_vector[i]
            box_second[i] = box_coord[i] + second_vector[i]
            box_third[i] = box_coord[i] + third_vector[i]
        surface[0].append(
            first_vector[0] * box_first[0] + first_vector[1] * box_first[1] + first_vector[2] * box_first[2])
        surface[1].append(
            first_vector[0] * box_coord[0] + first_vector[1] * box_coord[1] + first_vector[2] * box_coord[2])
        surface[2].append(
            second_vector[0] * box_second[0] + second_vector[1] * box_second[1] + second_vector[2] * box_second[2])
        surface[3].append(
            second_vector[0] * box_coord[0] + second_vector[1] * box_coord[1] + second_vector[2] * box_coord[2])
        surface[4].append(
            third_vector[0] * box_third[0] + third_vector[1] * box_third[1] + third_vector[2] * box_third[2])
        surface[5].append(
            third_vector[0] * box_coord[0] + third_vector[1] * box_coord[1] + third_vector[2] * box_coord[2])
        for i in range(6):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{self.surf_ids[1]}&{-self.surf_ids[2]}&{self.surf_ids[3]}&{-self.surf_ids[4]}' \
                   f'&{self.surf_ids[5]})'
        else:
            return f'({self.surf_ids[0]}:{-self.surf_ids[1]}:{self.surf_ids[2]}:{-self.surf_ids[3]}:{self.surf_ids[4]}' \
                   f':{-self.surf_ids[5]})'


class SPH(MacroBody):
    # 球体

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        if len(self.params) == 1:
            self.params = [0.0, 0.0, 0.0] + self.params
        try:
            surface = []
            surface_model = []
            self.surf_ids.append(max_surf_id + 1)
            surf_type = 'S'
            for i in range(0, 4):
                surface.append(self.params[i])
            surf = Surface.externalization(number=self.surf_ids[0], type=surf_type, parameters=surface,
                                        rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        except:
            print("SPH transfer error"+self.number)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'{-self.surf_ids[0]}'
        else:
            return f'{self.surf_ids[0]}'


class TORUS(MacroBody):
    # 圆环体

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface = []
        surface_model = []
        self.surf_ids.append(max_surf_id + 1)
        direction_vector = self.params[0:3]  # 方向向量
        center_coordinate = self.params[3:6]  # 中心点坐标
        circle_radius = self.params[6]  # 圆环半径
        tangent_radius = self.params[7]  # 切面半径
        if [direction_vector[0], direction_vector[1]] == [0, 0]:
            surf_type = 'TZ'
        elif [direction_vector[0], direction_vector[2]] == [0, 0]:
            surf_type = 'TY'
        elif [direction_vector[1], direction_vector[2]] == [0, 0]:
            surf_type = 'TX'
        for i in range(3):
            surface.append(center_coordinate[i])
        if sum(direction_vector) == -1:
            surface.append(-circle_radius + tangent_radius)
        else:
            surface.append(circle_radius - tangent_radius)
        surface.append(tangent_radius)
        surface.append(tangent_radius)

        surf = Surface.externalization(number=self.surf_ids[0], type=surf_type, parameters=surface,
                                       rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
        surf = surf.transfer()
        surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'{-self.surf_ids[0]}'
        else:
            return f'{self.surf_ids[0]}'


class ELL(MacroBody):
    # 旋转椭球面

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface = []
        surface_model = []
        self.surf_ids.append(max_surf_id + 1)
        if self.params[6] < 0:
            co_bottom = [self.params[0], self.params[1], self.params[2]]
            long_vector = [self.params[3], self.params[4], self.params[5]]
            short_vector_length = abs(self.params[6])
            long_vector_length = np.linalg.norm(long_vector)
            long_vector = [self.params[3] / long_vector_length, self.params[4] / long_vector_length, self.params[5] /
                           long_vector_length]
            a = long_vector_length * long_vector_length
            b = short_vector_length * short_vector_length
        else:
            a = self.params[6] * self.params[6]
            co_bottom = [(self.params[3] + self.params[0]) / 2, (self.params[4] + self.params[1]) / 2,
                         (self.params[5] + self.params[2]) / 2]
            long_vector = [self.params[3] - self.params[0], self.params[4] - self.params[1],
                           self.params[5] - self.params[2]]
            c = np.linalg.norm(long_vector) / 2
            long_vector = [long_vector[0] / (2 * c), long_vector[1] / (2 * c), long_vector[2] / (2 * c)]
            b = a - c * c
        c = -a + b
        du = a + long_vector[0] * long_vector[0] * c
        dv = a + long_vector[1] * long_vector[1] * c
        dw = a + long_vector[2] * long_vector[2] * c
        # return SQ
        if (long_vector[0] == 0 or long_vector[1] == 0) and (long_vector[0] == 0 or long_vector[2] == 0) and (
                long_vector[1] == 0 or long_vector[2] == 0):
            surface.append(du / (a * b))
            surface.append(dv / (a * b))
            surface.append(dw / (a * b))
            surface.append(0)
            surface.append(0)
            surface.append(0)
            surface.append(-1)
            surface.append(co_bottom[0])
            surface.append(co_bottom[1])
            surface.append(co_bottom[2])
            surf_type = 'SQ'
        # return GQ
        else:
            surface.append(du / (a * b))
            surface.append(dv / (a * b))
            surface.append(dw / (a * b))
            surface.append((2 * long_vector[0] * long_vector[1] * c) / (a * b))
            surface.append((2 * long_vector[1] * long_vector[2] * c) / (a * b))
            surface.append((2 * long_vector[0] * long_vector[2] * c) / (a * b))
            surface.append((-2 * co_bottom[0] * du - 2 * long_vector[0] * long_vector[1] * co_bottom[1] * c - 2 *
                            long_vector[0] * long_vector[2] * co_bottom[2] * c) / (a * b))
            surface.append((-2 * co_bottom[1] * dv - 2 * long_vector[0] * long_vector[1] * co_bottom[0] * c - 2 *
                            long_vector[1] * long_vector[2] * co_bottom[2] * c) / (a * b))
            surface.append((-2 * co_bottom[2] * dw - 2 * long_vector[0] * long_vector[2] * co_bottom[0] * c - 2 *
                            long_vector[1] * long_vector[2] * co_bottom[1] * c) / (a * b))
            surface.append((co_bottom[0] * co_bottom[0] * du + co_bottom[1] * co_bottom[1] * dv + co_bottom[2] *
                            co_bottom[2] * dw + co_bottom[0] * co_bottom[1] * surface[3] * a * b + co_bottom[1] *
                            co_bottom[2] * surface[4] * a * b - a * b + co_bottom[0] * co_bottom[2] * surface[
                                5] * a * b)
                           / (a * b))
            surf_type = 'GQ'
        surf = Surface.externalization(number=self.surf_ids[0], type=surf_type, parameters=surface,
                                       rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
        surf = surf.transfer()
        surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'{-self.surf_ids[0]}'
        else:
            return f'{self.surf_ids[0]}'

class HEX(MacroBody):
    # 六棱柱

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface_model = []
        surface = [[] for _ in range(8)]
        co_bottom = [self.params[0], self.params[1], self.params[2]]  # 获取底面坐标
        axial_vector = [self.params[3], self.params[4], self.params[5]]  # 轴向向量
        axial_length = np.linalg.norm(axial_vector)  # 圆柱轴向长度
        axial_vector = [axial_vector[0] / axial_length, axial_vector[1] / axial_length,
                        axial_vector[2] / axial_length]  # 轴向单位向量
        vector_r = [self.params[6], self.params[7], self.params[8]]  # 中心到第一个面的向量
        vector_s = np.zeros(3)
        vector_t = np.zeros(3)
        surf_type = ['P', 'P', 'P', 'P', 'P', 'P', 'P', 'P']
        if len(self.params) == 9:
            # 这里参数长度为9的情况是正六棱柱
            vector_s = [
                vector_r[0] / 2 + (math.sqrt(3) / 2) * (vector_r[1] * axial_vector[2] - vector_r[2] * axial_vector[1]),
                vector_r[1] / 2 + (math.sqrt(3) / 2) * (vector_r[2] * axial_vector[0] - vector_r[0] * axial_vector[2]),
                vector_r[2] / 2 + (math.sqrt(3) / 2) * (vector_r[0] * axial_vector[1] - vector_r[1] * axial_vector[0])]
            vector_t = [
                vector_r[0] / 2 - (math.sqrt(3) / 2) * (vector_r[1] * axial_vector[2] - vector_r[2] * axial_vector[1]),
                vector_r[1] / 2 - (math.sqrt(3) / 2) * (vector_r[2] * axial_vector[0] - vector_r[0] * axial_vector[2]),
                vector_r[2] / 2 - (math.sqrt(3) / 2) * (vector_r[0] * axial_vector[1] - vector_r[1] * axial_vector[0])]
        else:
            vector_s = [self.params[9], self.params[10], self.params[11]]
            vector_t = [self.params[12], self.params[13], self.params[14]]
        for i in range(3):
            surface[0].append(vector_r[i])
            surface[1].append(vector_r[i])
            surface[2].append(vector_s[i])
            surface[3].append(vector_s[i])
            surface[4].append(vector_t[i])
            surface[5].append(vector_t[i])
            surface[6].append(axial_vector[i])
            surface[7].append(axial_vector[i])
        surface[0].append(vector_r[0] * (co_bottom[0] + vector_r[0]) + vector_r[1] * (co_bottom[1] + vector_r[1]) +
                          vector_r[2] * (co_bottom[2] + vector_r[2]))
        surface[1].append(vector_r[0] * (co_bottom[0] - vector_r[0]) + vector_r[1] * (co_bottom[1] - vector_r[1]) +
                          vector_r[2] * (co_bottom[2] - vector_r[2]))
        surface[2].append(vector_s[0] * (co_bottom[0] + vector_s[0]) + vector_s[1] * (co_bottom[1] + vector_s[1]) +
                          vector_s[2] * (co_bottom[2] + vector_s[2]))
        surface[3].append(vector_s[0] * (co_bottom[0] - vector_s[0]) + vector_s[1] * (co_bottom[1] - vector_s[1]) +
                          vector_s[2] * (co_bottom[2] - vector_s[2]))
        surface[4].append(vector_t[0] * (co_bottom[0] + vector_t[0]) + vector_t[1] * (co_bottom[1] + vector_t[1]) +
                          vector_t[2] * (co_bottom[2] + vector_t[2]))
        surface[5].append(vector_t[0] * (co_bottom[0] - vector_t[0]) + vector_t[1] * (co_bottom[1] - vector_t[1]) +
                          vector_t[2] * (co_bottom[2] - vector_t[2]))
        surface[6].append(axial_vector[0] * (co_bottom[0] + axial_vector[0] * axial_length) + axial_vector[1] *
                          (co_bottom[1] + axial_vector[1] * axial_length) + axial_vector[2] * (
                                  co_bottom[2] + axial_vector[2] *
                                  axial_length))
        surface[7].append(
            axial_vector[0] * co_bottom[0] + axial_vector[1] * co_bottom[1] + axial_vector[2] * co_bottom[2])
        for i in range(0, 8):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{self.surf_ids[1]}&{-self.surf_ids[2]}&{self.surf_ids[3]}&{-self.surf_ids[4]}' \
                   f'&{self.surf_ids[5]}&{-self.surf_ids[6]}&{self.surf_ids[7]})'
        else:
            return f'({self.surf_ids[0]}:{-self.surf_ids[1]}:{self.surf_ids[2]}:{-self.surf_ids[3]}:{self.surf_ids[4]}' \
                   f':{-self.surf_ids[5]}:{self.surf_ids[6]}:{-self.surf_ids[7]})'
                   

class HEX(MacroBody):
    # 六棱柱

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface_model = []
        surface = [[] for _ in range(8)]
        co_bottom = [self.params[0], self.params[1], self.params[2]]  # 获取底面坐标
        axial_vector = [self.params[3], self.params[4], self.params[5]]  # 轴向向量
        axial_length = np.linalg.norm(axial_vector)  # 圆柱轴向长度
        axial_vector = [axial_vector[0] / axial_length, axial_vector[1] / axial_length,
                        axial_vector[2] / axial_length]  # 轴向单位向量
        vector_r = [self.params[6], self.params[7], self.params[8]]  # 中心到第一个面的向量
        vector_s = np.zeros(3)
        vector_t = np.zeros(3)
        surf_type = ['P', 'P', 'P', 'P', 'P', 'P', 'P', 'P']
        if len(self.params) == 9:
            # 这里参数长度为9的情况是正六棱柱
            vector_s = [
                vector_r[0] / 2 + (math.sqrt(3) / 2) * (vector_r[1] * axial_vector[2] - vector_r[2] * axial_vector[1]),
                vector_r[1] / 2 + (math.sqrt(3) / 2) * (vector_r[2] * axial_vector[0] - vector_r[0] * axial_vector[2]),
                vector_r[2] / 2 + (math.sqrt(3) / 2) * (vector_r[0] * axial_vector[1] - vector_r[1] * axial_vector[0])]
            vector_t = [
                vector_r[0] / 2 - (math.sqrt(3) / 2) * (vector_r[1] * axial_vector[2] - vector_r[2] * axial_vector[1]),
                vector_r[1] / 2 - (math.sqrt(3) / 2) * (vector_r[2] * axial_vector[0] - vector_r[0] * axial_vector[2]),
                vector_r[2] / 2 - (math.sqrt(3) / 2) * (vector_r[0] * axial_vector[1] - vector_r[1] * axial_vector[0])]
        else:
            vector_s = [self.params[9], self.params[10], self.params[11]]
            vector_t = [self.params[12], self.params[13], self.params[14]]
        for i in range(3):
            surface[0].append(vector_r[i])
            surface[1].append(vector_r[i])
            surface[2].append(vector_s[i])
            surface[3].append(vector_s[i])
            surface[4].append(vector_t[i])
            surface[5].append(vector_t[i])
            surface[6].append(axial_vector[i])
            surface[7].append(axial_vector[i])
        surface[0].append(vector_r[0] * (co_bottom[0] + vector_r[0]) + vector_r[1] * (co_bottom[1] + vector_r[1]) +
                          vector_r[2] * (co_bottom[2] + vector_r[2]))
        surface[1].append(vector_r[0] * (co_bottom[0] - vector_r[0]) + vector_r[1] * (co_bottom[1] - vector_r[1]) +
                          vector_r[2] * (co_bottom[2] - vector_r[2]))
        surface[2].append(vector_s[0] * (co_bottom[0] + vector_s[0]) + vector_s[1] * (co_bottom[1] + vector_s[1]) +
                          vector_s[2] * (co_bottom[2] + vector_s[2]))
        surface[3].append(vector_s[0] * (co_bottom[0] - vector_s[0]) + vector_s[1] * (co_bottom[1] - vector_s[1]) +
                          vector_s[2] * (co_bottom[2] - vector_s[2]))
        surface[4].append(vector_t[0] * (co_bottom[0] + vector_t[0]) + vector_t[1] * (co_bottom[1] + vector_t[1]) +
                          vector_t[2] * (co_bottom[2] + vector_t[2]))
        surface[5].append(vector_t[0] * (co_bottom[0] - vector_t[0]) + vector_t[1] * (co_bottom[1] - vector_t[1]) +
                          vector_t[2] * (co_bottom[2] - vector_t[2]))
        surface[6].append(axial_vector[0] * (co_bottom[0] + axial_vector[0] * axial_length) + axial_vector[1] *
                          (co_bottom[1] + axial_vector[1] * axial_length) + axial_vector[2] * (
                                  co_bottom[2] + axial_vector[2] *
                                  axial_length))
        surface[7].append(
            axial_vector[0] * co_bottom[0] + axial_vector[1] * co_bottom[1] + axial_vector[2] * co_bottom[2])
        for i in range(0, 8):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{self.surf_ids[1]}&{-self.surf_ids[2]}&{self.surf_ids[3]}&{-self.surf_ids[4]}' \
                   f'&{self.surf_ids[5]}&{-self.surf_ids[6]}&{self.surf_ids[7]})'
        else:
            return f'({self.surf_ids[0]}:{-self.surf_ids[1]}:{self.surf_ids[2]}:{-self.surf_ids[3]}:{self.surf_ids[4]}' \
                   f':{-self.surf_ids[5]}:{self.surf_ids[6]}:{-self.surf_ids[7]})'


class REC(MacroBody):
    # 椭圆柱

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surf_type = []
        surface_model = []
        surface = [[] for _ in range(3)]
        sho_len = 0
        # 这里参数可以为10个或者12个
        if len(self.params) == 10:
            short_vec = self.params[9]
            sho_len = abs(short_vec)
        elif len(self.params) == 12:
            short_vec = self.params[9:12]
            sho_len = np.linalg.norm(short_vec)
        co_bottom = self.params[0:3]  # 底面坐标
        axial_vec = self.params[3:6]  # 轴向向量
        ax_len = np.linalg.norm(axial_vec)  # 轴向长度
        long_vec = self.params[6:9]
        lo_len = np.linalg.norm(long_vec)  # 长轴长度
        a2 = lo_len * lo_len
        b2 = sho_len * sho_len
        h2 = ax_len * ax_len
        # return SQ
        if (axial_vec[0] == 0 or axial_vec[1] == 0) and (axial_vec[0] == 0 or axial_vec[2] == 0) and (
                axial_vec[1] == 0 or axial_vec[2] == 0):
            surface[0].append((a2 * (1 - axial_vec[0] * axial_vec[0] / h2) - long_vec[0] * long_vec[0] * (1 - b2 / a2))
                              / (a2 * b2))
            surface[0].append((a2 * (1 - axial_vec[1] * axial_vec[1] / h2) - long_vec[1] * long_vec[1] * (1 - b2 / a2))
                              / (a2 * b2))
            surface[0].append((a2 * (1 - axial_vec[2] * axial_vec[2] / h2) - long_vec[2] * long_vec[2] * (1 - b2 / a2))
                              / (a2 * b2))
            surface[0].append(0)
            surface[0].append(0)
            surface[0].append(0)
            if [axial_vec[0], axial_vec[1]] == [0, 0]:
                surf_type = ['SQ', 'P', 'P']
                surface[0][2] = 0
                r = co_bottom[0] * co_bottom[0] + co_bottom[1] * co_bottom[1]
                surface[0].append((a2 * r - (co_bottom[0] * co_bottom[0] * long_vec[0] * long_vec[0] + co_bottom[1] *
                                             co_bottom[1] * long_vec[1] * long_vec[1]) * (1 - b2 / a2) - a2 * b2 -
                                   surface[0][0] * a2 * b2 * co_bottom[0] * co_bottom[0] - surface[0][1] * a2 * b2 *
                                   co_bottom[1]
                                   * co_bottom[1]) / (a2 * b2))
            elif [axial_vec[0], axial_vec[2]] == [0, 0]:
                surf_type = ['SQ', 'P', 'P']
                surface[0][1] = 0
                r = co_bottom[0] * co_bottom[0] + co_bottom[2] * co_bottom[2]
                surface[0].append((a2 * r - (co_bottom[0] * co_bottom[0] * long_vec[0] * long_vec[0] + co_bottom[2] *
                                             co_bottom[2] * long_vec[2] * long_vec[2]) * (1 - b2 / a2) - a2 * b2 -
                                   surface[0][0] * a2 * b2 * co_bottom[0] * co_bottom[0] - surface[0][2] * a2 * b2 *
                                   co_bottom[2]
                                   * co_bottom[2]) / (a2 * b2))
            elif [axial_vec[1], axial_vec[2]] == [0, 0]:
                surf_type = ['SQ', 'P', 'P']
                surface[0][0] = 0
                r = co_bottom[1] * co_bottom[1] + co_bottom[2] * co_bottom[2]
                surface[0].append((a2 * r - (co_bottom[1] * co_bottom[1] * long_vec[1] * long_vec[1] + co_bottom[2] *
                                             co_bottom[2] * long_vec[2] * long_vec[2]) * (1 - b2 / a2) - a2 * b2 -
                                   surface[0][1] * a2 * b2 * co_bottom[1] * co_bottom[1] - surface[0][2] * a2 * b2 *
                                   co_bottom[2]
                                   * co_bottom[2]) / (a2 * b2))
            surface[0].append(co_bottom[0])
            surface[0].append(co_bottom[1])
            surface[0].append(co_bottom[2])

        # return GQ
        else:
            surf_type = ['GQ', 'P', 'P']
            surface[0].append((a2 * (1 - axial_vec[0] * axial_vec[0] / h2) - long_vec[0] * long_vec[0] * (1 - b2 / a2))
                              / (a2 * b2))
            surface[0].append((a2 * (1 - axial_vec[1] * axial_vec[1] / h2) - long_vec[1] * long_vec[1] * (1 - b2 / a2))
                              / (a2 * b2))
            surface[0].append((a2 * (1 - axial_vec[2] * axial_vec[2] / h2) - long_vec[2] * long_vec[2] * (1 - b2 / a2))
                              / (a2 * b2))
            surface[0].append(
                (-2 * a2 * axial_vec[0] * axial_vec[1] / h2 - 2 * long_vec[0] * long_vec[1] * (1 - b2 / a2))
                / (a2 * b2))
            surface[0].append(
                (-2 * a2 * axial_vec[1] * axial_vec[2] / h2 - 2 * long_vec[1] * long_vec[2] * (1 - b2 / a2))
                / (a2 * b2))
            surface[0].append(
                (-2 * a2 * axial_vec[0] * axial_vec[2] / h2 - 2 * long_vec[0] * long_vec[2] * (1 - b2 / a2))
                / (a2 * b2))
            surface[0].append(
                (2 * a2 * ((axial_vec[0] / h2) * (co_bottom[1] * axial_vec[1] + co_bottom[2] * axial_vec[2])
                           - co_bottom[0] * ((1 - axial_vec[0] * axial_vec[0] / h2) -
                                             (long_vec[0] * long_vec[0] / a2) * (1 - b2 / a2))) + 2 * (1 - b2 / a2)
                 * (co_bottom[1] * long_vec[0] * long_vec[1] + co_bottom[2] * long_vec[0] * long_vec[2]))
                / (a2 * b2))
            surface[0].append(
                (2 * a2 * ((axial_vec[1] / h2) * (co_bottom[0] * axial_vec[0] + co_bottom[2] * axial_vec[2])
                           - co_bottom[1] * ((1 - axial_vec[1] * axial_vec[1] / h2) -
                                             (long_vec[1] * long_vec[1] / a2) * (1 - b2 / a2))) + 2 * (1 - b2 / a2)
                 * (co_bottom[0] * long_vec[0] * long_vec[1] + co_bottom[2] * long_vec[1] * long_vec[2]))
                / (a2 * b2))
            surface[0].append(
                (2 * a2 * ((axial_vec[2] / h2) * (co_bottom[0] * axial_vec[0] + co_bottom[1] * axial_vec[1])
                           - co_bottom[2] * ((1 - axial_vec[2] * axial_vec[2] / h2) -
                                             (long_vec[2] * long_vec[2] / a2) * (1 - b2 / a2))) + 2 * (1 - b2 / a2)
                 * (co_bottom[0] * long_vec[0] * long_vec[2] + co_bottom[1] * long_vec[1] * long_vec[2]))
                / (a2 * b2))
            surface[0].append((a2 * ((-2 * (co_bottom[0] * co_bottom[1] * axial_vec[0] * axial_vec[1] + co_bottom[1] *
                                            co_bottom[2] * axial_vec[1] * axial_vec[2] + co_bottom[0] * co_bottom[2] *
                                            axial_vec[0] * axial_vec[2]) / h2) + co_bottom[0] * co_bottom[0] *
                                     (1 - axial_vec[0] * axial_vec[0] / h2) + co_bottom[1] * co_bottom[1] *
                                     (1 - axial_vec[1] * axial_vec[1] / h2) + co_bottom[2] * co_bottom[2] *
                                     (1 - axial_vec[2] * axial_vec[2] / h2)) -
                               (co_bottom[0] * co_bottom[0] * long_vec[0] * long_vec[0] +
                                co_bottom[1] * co_bottom[1] * long_vec[1] * long_vec[1] +
                                co_bottom[2] * co_bottom[2] * long_vec[2] * long_vec[2]) * (1 - b2 / a2) - a2 * b2 -
                               2 * (co_bottom[0] * co_bottom[1] * long_vec[0] * long_vec[1] + co_bottom[1] *
                                    co_bottom[2] * long_vec[1] * long_vec[2] + co_bottom[0] * co_bottom[2] * long_vec[0]
                                    * long_vec[2]) * (1 - b2 / a2)) / (a2 * b2))
        surface[1].append(axial_vec[0])
        surface[1].append(axial_vec[1])
        surface[1].append(axial_vec[2])
        surface[1].append(axial_vec[0] * (co_bottom[0] + axial_vec[0]) + axial_vec[1] * (co_bottom[1] + axial_vec[1]) +
                          axial_vec[2] * (co_bottom[2] + axial_vec[2]))
        surface[2].append(axial_vec[0])
        surface[2].append(axial_vec[1])
        surface[2].append(axial_vec[2])
        surface[2].append(axial_vec[0] * co_bottom[0] + axial_vec[1] * co_bottom[1] + axial_vec[2] * co_bottom[2])
        for i in range(0, 3):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{-self.surf_ids[1]}&{self.surf_ids[2]})'
        else:
            return f'({self.surf_ids[0]}:{self.surf_ids[1]}:{-self.surf_ids[2]})'


class TRC(MacroBody):
    # 圆台

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surf_type = []
        surface_parameter = [[] for _ in range(3)]
        surface_model = []
        co_bottom = [self.params[0], self.params[1], self.params[2]]
        axial_vector = [self.params[3], self.params[4], self.params[5]]
        axial_vector_length = np.linalg.norm(axial_vector)
        axial_vector = [axial_vector[0] / axial_vector_length, axial_vector[1] / axial_vector_length, axial_vector[2] /
                        axial_vector_length]
        out_radius = self.params[6]  # 大圆半径
        inner_radius = self.params[7]  # 小圆半径
        # 轴平行于坐标轴
        if (axial_vector[0] == 0 or axial_vector[1] == 0) and (axial_vector[0] == 0 or axial_vector[2] == 0) and (
                axial_vector[1] == 0 or axial_vector[2] == 0):
            surface = np.zeros(5)
            surface[3] = math.pow((out_radius - inner_radius) / axial_vector_length, 2)
            if [axial_vector[0], axial_vector[1]] == [0, 0]:
                surf_type = ['K/Z', 'P', 'P']
                if axial_vector[2] > 0:
                    surface[4] = -1
                    surface[2] = co_bottom[2] + out_radius / math.sqrt(surface[3])
                elif axial_vector[2] < 0:
                    surface[4] = 1
                    surface[2] = co_bottom[2] - out_radius / math.sqrt(surface[3])
                surface[0] = co_bottom[0]
                surface[1] = co_bottom[1]
            elif [axial_vector[0], axial_vector[2]] == [0, 0]:
                surf_type = ['K/Y', 'P', 'P']
                if axial_vector[1] > 0:
                    surface[4] = -1
                    surface[1] = co_bottom[1] + out_radius / math.sqrt(surface[3])
                elif axial_vector[1] < 0:
                    surface[4] = 1
                    surface[1] = co_bottom[1] - out_radius / math.sqrt(surface[3])
                surface[0] = co_bottom[0]
                surface[2] = co_bottom[2]
            elif [axial_vector[1], axial_vector[2]] == [0, 0]:
                surf_type = ['K/X', 'P', 'P']
                if axial_vector[0] > 0:
                    surface[4] = -1
                    surface[0] = co_bottom[0] + out_radius / math.sqrt(surface[3])
                elif axial_vector[0] < 0:
                    surface[4] = 1
                    surface[0] = co_bottom[0] - out_radius / math.sqrt(surface[3])
                surface[1] = co_bottom[1]
                surface[2] = co_bottom[2]
            for i in range(5):
                surface_parameter[0].append(surface[i])
        # 轴不平行坐标轴
        else:
            surf_type = ['GQ', 'P', 'P']
            surface = np.zeros(10)
            # return gq
            surface[0] = 1 - axial_vector[0] * axial_vector[0]
            surface[1] = 1 - axial_vector[1] * axial_vector[1]
            surface[2] = 1 - axial_vector[2] * axial_vector[2]
            surface[3] = -2 * axial_vector[0] * axial_vector[1]
            surface[4] = -2 * axial_vector[1] * axial_vector[2]
            surface[5] = -2 * axial_vector[0] * axial_vector[2]
            surface[6] = -co_bottom[1] * surface[3] - co_bottom[2] * surface[5] - 2 * co_bottom[0] * surface[0]
            surface[7] = -co_bottom[0] * surface[3] - co_bottom[2] * surface[4] - 2 * co_bottom[1] * surface[1]
            surface[8] = -co_bottom[0] * surface[5] - co_bottom[1] * surface[4] - 2 * co_bottom[2] * surface[2]
            surface[9] = co_bottom[0] * co_bottom[1] * surface[3] + co_bottom[1] * co_bottom[2] * surface[4] + \
                         co_bottom[0] * co_bottom[2] * surface[5] + co_bottom[0] * co_bottom[0] * surface[0] + \
                         co_bottom[1] * co_bottom[1] * surface[1] + co_bottom[2] * co_bottom[2] * surface[2] - \
                         out_radius * out_radius
            # add terms to change cylinder to cone
            t2 = (out_radius - inner_radius) * (out_radius - inner_radius) / (axial_vector_length * axial_vector_length)
            dp = -(out_radius / inner_radius) * axial_vector_length / ((out_radius / inner_radius) - 1)
            surface[0] = surface[0] - t2 * axial_vector[0] * axial_vector[0]
            surface[1] = surface[1] - t2 * axial_vector[1] * axial_vector[1]
            surface[2] = surface[2] - t2 * axial_vector[2] * axial_vector[2]
            surface[6] = surface[6] - t2 * (
                    co_bottom[1] * surface[3] + co_bottom[2] * surface[5] - 2 * co_bottom[0] * axial_vector[0]
                    * axial_vector[0] + 2 * dp * axial_vector[0])
            surface[7] = surface[7] - t2 * (
                    co_bottom[0] * surface[3] + co_bottom[2] * surface[4] - 2 * co_bottom[1] * axial_vector[1]
                    * axial_vector[1] + 2 * dp * axial_vector[1])
            surface[8] = surface[8] - t2 * (
                    co_bottom[0] * surface[5] + co_bottom[1] * surface[4] - 2 * co_bottom[2] * axial_vector[2]
                    * axial_vector[2] + 2 * dp * axial_vector[2])
            surface[9] = surface[9] + out_radius * out_radius - t2 * (-co_bottom[0] * co_bottom[1] * surface[3] -
                                                                      co_bottom[1] * co_bottom[2] * surface[4] -
                                                                      co_bottom[0] * co_bottom[2] * surface[5] +
                                                                      co_bottom[0] * co_bottom[0] * axial_vector[0] *
                                                                      axial_vector[0] + co_bottom[1] * co_bottom[1] *
                                                                      axial_vector[1] * axial_vector[1] + co_bottom[2] *
                                                                      co_bottom[2] * axial_vector[2] * axial_vector[2] +
                                                                      dp * dp - 2 * dp * (co_bottom[0] * axial_vector[0]
                                                                                          + co_bottom[1] * axial_vector[
                                                                                              1] +
                                                                                          co_bottom[2] * axial_vector[
                                                                                              2]))
            surface[3] = surface[3] - 2 * t2 * axial_vector[0] * axial_vector[1]
            surface[4] = surface[4] - 2 * t2 * axial_vector[1] * axial_vector[2]
            surface[5] = surface[5] - 2 * t2 * axial_vector[0] * axial_vector[2]
            for i in range(10):
                surface_parameter[0].append(surface[i])
        surface_parameter[1].append(axial_vector[0])
        surface_parameter[1].append(axial_vector[1])
        surface_parameter[1].append(axial_vector[2])
        surface_parameter[1].append(
            axial_vector[0] * (co_bottom[0] + axial_vector[0] * axial_vector_length) + axial_vector[1] *
            (co_bottom[1] + axial_vector[1] * axial_vector_length) + axial_vector[2] *
            (co_bottom[2] + axial_vector[2] * axial_vector_length))
        surface_parameter[2].append(axial_vector[0])
        surface_parameter[2].append(axial_vector[1])
        surface_parameter[2].append(axial_vector[2])
        surface_parameter[2].append(axial_vector[0] * co_bottom[0] + axial_vector[1] * co_bottom[1] + axial_vector[2]
                                    * co_bottom[2])
        for i in range(0, 3):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface_parameter[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            return f'({-self.surf_ids[0]}&{-self.surf_ids[1]}&{self.surf_ids[2]})'
        else:
            return f'({self.surf_ids[0]}:{self.surf_ids[1]}:{-self.surf_ids[2]})'


class WED(MacroBody):
    # 三棱柱
    surface_sign = []  # 判断在该法向量下，侧面三个平面的符号


    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface_model = []
        surface = [[] for _ in range(5)]
        co_bottom = self.params[0:3]  # 底面一个顶点坐标
        axial_vector = self.params[9:12]  # 轴向向量
        first_side_vector = self.params[3:6]  # 底面第一个边向量
        second_side_vector = self.params[6:9]  # 底面第二个边向量
        first_normal_vector = np.zeros(3)  # 包含第一条底边侧面的法向量
        second_normal_vector = np.zeros(3)  # 包含第二条底边侧面的法向量
        third_normal_vector = np.zeros(3)  # 包含第三条底边侧面的法向量
        third_side_vector = np.zeros(3)  # 底面第三个边向量
        surf_type = ['P', 'P', 'P', 'P', 'P']
        for i in range(0, 3):
            third_side_vector[i] = first_side_vector[i] - second_side_vector[i]
        first_normal_vector[0] = axial_vector[1] * first_side_vector[2] - first_side_vector[1] * axial_vector[2]
        first_normal_vector[1] = axial_vector[2] * first_side_vector[0] - first_side_vector[2] * axial_vector[0]
        first_normal_vector[2] = -axial_vector[1] * first_side_vector[0] + first_side_vector[1] * axial_vector[0]
        second_normal_vector[0] = second_side_vector[1] * axial_vector[2] - axial_vector[1] * second_side_vector[2]
        second_normal_vector[1] = axial_vector[2] * second_side_vector[0] - second_side_vector[2] * axial_vector[0]
        second_normal_vector[2] = -axial_vector[1] * second_side_vector[0] + second_side_vector[1] * axial_vector[0]
        third_normal_vector[0] = axial_vector[1] * third_side_vector[2] - third_side_vector[1] * axial_vector[2]
        third_normal_vector[1] = axial_vector[2] * third_side_vector[0] - third_side_vector[2] * axial_vector[0]
        third_normal_vector[2] = -axial_vector[1] * third_side_vector[0] + third_side_vector[1] * axial_vector[0]
        for i in range(3):
            surface[0].append(third_normal_vector[i])
            surface[1].append(second_normal_vector[i])
            surface[2].append(first_normal_vector[i])
            surface[3].append(axial_vector[i])
            surface[4].append(axial_vector[i])
        surface[0].append(third_normal_vector[0] * (co_bottom[0] + first_side_vector[0])
                          + third_normal_vector[1] * (co_bottom[1] + first_side_vector[1])
                          + third_normal_vector[2] * (co_bottom[2] + first_side_vector[2]))
        surface[1].append(second_normal_vector[0] * co_bottom[0] + second_normal_vector[1] * co_bottom[1] +
                          second_normal_vector[2] * co_bottom[2])
        surface[2].append(first_normal_vector[0] * co_bottom[0] + first_normal_vector[1] * co_bottom[1] +
                          first_normal_vector[2] * co_bottom[2])
        surface[3].append(axial_vector[0] * (co_bottom[0] + axial_vector[0]) +
                          axial_vector[1] * (co_bottom[1] + axial_vector[1]) +
                          axial_vector[2] * (co_bottom[2] + axial_vector[2]))
        surface[4].append(
            axial_vector[0] * co_bottom[0] + axial_vector[1] * co_bottom[1] + axial_vector[2] * co_bottom[2])
        self.surface_sign.append(
            -third_normal_vector[0] * first_side_vector[0] - third_normal_vector[1] * first_side_vector[1]
            - third_normal_vector[2] * first_side_vector[2])
        self.surface_sign.append(
            second_normal_vector[0] * first_side_vector[0] + second_normal_vector[1] * first_side_vector[1]
            + second_normal_vector[2] * first_side_vector[2])
        self.surface_sign.append(
            first_normal_vector[0] * second_side_vector[0] + first_normal_vector[1] * second_side_vector[1]
            + first_normal_vector[2] * second_side_vector[2])
        for i in range(0, 5):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            bool_list = '('
            for i, value in enumerate(self.surface_sign):
                if value < 0:
                    bool_list += f'{-self.surf_ids[i]}&'
                else:
                    bool_list += f'{self.surf_ids[i]}&'
            bool_list += f'{-self.surf_ids[3]}&{self.surf_ids[4]})'
            return bool_list
        else:
            bool_list = '('
            for i, value in enumerate(self.surface_sign):
                if value < 0:
                    bool_list += f'{self.surf_ids[i]}:'
                else:
                    bool_list += f'{-self.surf_ids[i]}:'
            bool_list += f'{self.surf_ids[3]}:{-self.surf_ids[4]})'
            return bool_list


class SEC(MacroBody):
    # 圆柱扇体
    surface_sign = []  # 判断在该法向量下，圆柱扇体两个侧平面内部与外部的符号

    def __init__(self, number, type, params, rotate=None, rotate_angle=None, move=None):
        super().__init__(number=number, type=type, params=params, rotate=rotate, rotate_angle=rotate_angle, move=move)

    def transfer(self, max_surf_id=0):
        surface_model = []
        surf_type = []
        surface = [[] for _ in range(6)]
        co_bottom = self.params[0:3]  # 底面坐标
        axial_vector = self.params[3:6]  # 轴向向量
        inner_radius = self.params[6]  # 内径
        out_radius = self.params[7]  # 外径
        inner_theta = self.params[8]  # 内侧面与坐标轴夹角
        out_theta = self.params[9]  # 外侧面与坐标轴夹角
        med_theta = (inner_theta + out_theta) / 2
        med_radius = (inner_radius + out_radius) / 2
        in_coordinate = np.zeros(3)  # 圆柱扇体中一点的坐标
        vector_p1 = np.zeros(3)
        vector_p2 = np.zeros(3)
        first_normal_vector = np.zeros(3)
        second_normal_vector = np.zeros(3)
        # return p
        surface[0].append(axial_vector[0])
        surface[0].append(axial_vector[1])
        surface[0].append(axial_vector[2])
        surface[0].append(
            axial_vector[0] * co_bottom[0] + axial_vector[1] * co_bottom[1] + axial_vector[2] * co_bottom[2])
        # return p
        surface[1].append(axial_vector[0])
        surface[1].append(axial_vector[1])
        surface[1].append(axial_vector[2])
        surface[1].append(
            axial_vector[0] * (co_bottom[0] + axial_vector[0]) + axial_vector[1] * (co_bottom[1] + axial_vector[1])
            + axial_vector[2] * (co_bottom[2] + axial_vector[2]))
        if [axial_vector[0], axial_vector[1]] == [0, 0]:
            surf_type = ['P', 'P', 'C/Z', 'C/Z', 'P', 'P']
            surface[2].append(co_bottom[0])
            surface[2].append(co_bottom[1])
            surface[2].append(inner_radius)
            surface[3].append(co_bottom[0])
            surface[3].append(co_bottom[1])
            surface[3].append(out_radius)
            in_coordinate[0] = med_radius * math.cos(med_theta * 2 * math.pi / 360)
            in_coordinate[1] = med_radius * math.sin(med_theta * 2 * math.pi / 360)
            in_coordinate[2] = 0
            vector_p1[0] = math.cos(inner_theta * 2 * math.pi / 360)
            vector_p1[1] = math.sin(inner_theta * 2 * math.pi / 360)
            vector_p1[2] = 0
            vector_p2[0] = math.cos(out_theta * 2 * math.pi / 360)
            vector_p2[1] = math.sin(out_theta * 2 * math.pi / 360)
            vector_p2[2] = 0
        elif [axial_vector[0], axial_vector[2]] == [0, 0]:
            surf_type = ['P', 'P', 'C/Y', 'C/Y', 'P', 'P']
            surface[2].append(co_bottom[0])
            surface[2].append(co_bottom[2])
            surface[2].append(inner_radius)
            surface[3].append(co_bottom[0])
            surface[3].append(co_bottom[2])
            surface[3].append(out_radius)
            in_coordinate[2] = med_radius * math.cos(med_theta * 2 * math.pi / 360)
            in_coordinate[0] = med_radius * math.sin(med_theta * 2 * math.pi / 360)
            in_coordinate[1] = 0
            vector_p1[2] = math.cos(inner_theta * 2 * math.pi / 360)
            vector_p1[0] = math.sin(inner_theta * 2 * math.pi / 360)
            vector_p1[1] = 0
            vector_p2[2] = math.cos(out_theta * 2 * math.pi / 360)
            vector_p2[0] = math.sin(out_theta * 2 * math.pi / 360)
            vector_p2[1] = 0
        elif [axial_vector[1], axial_vector[2]] == [0, 0]:
            surf_type = ['P', 'P', 'C/X', 'C/X', 'P', 'P']
            surface[2].append(co_bottom[1])
            surface[2].append(co_bottom[2])
            surface[2].append(inner_radius)
            surface[3].append(co_bottom[1])
            surface[3].append(co_bottom[2])
            surface[3].append(out_radius)
            in_coordinate[1] = med_radius * math.cos(med_theta * 2 * math.pi / 360)
            in_coordinate[2] = med_radius * math.sin(med_theta * 2 * math.pi / 360)
            in_coordinate[0] = 0
            vector_p1[1] = math.cos(inner_theta * 2 * math.pi / 360)
            vector_p1[2] = math.sin(inner_theta * 2 * math.pi / 360)
            vector_p1[0] = 0
            vector_p2[1] = math.cos(out_theta * 2 * math.pi / 360)
            vector_p2[2] = math.sin(out_theta * 2 * math.pi / 360)
            vector_p2[0] = 0
        first_normal_vector[0] = axial_vector[1] * vector_p1[2] - vector_p1[1] * axial_vector[2]
        first_normal_vector[1] = axial_vector[2] * vector_p1[0] - vector_p1[2] * axial_vector[0]
        first_normal_vector[2] = axial_vector[0] * vector_p1[1] - vector_p1[0] * axial_vector[1]
        second_normal_vector[0] = vector_p2[2] * axial_vector[1] - axial_vector[2] * vector_p2[1]
        second_normal_vector[1] = axial_vector[2] * vector_p2[0] - vector_p2[2] * axial_vector[0]
        second_normal_vector[2] = axial_vector[0] * vector_p2[1] - vector_p2[0] * axial_vector[1]
        for i in range(3):
            surface[4].append(first_normal_vector[i])
            surface[5].append(second_normal_vector[i])
        self.surface_sign.append(first_normal_vector[0] * in_coordinate[0] + first_normal_vector[1] * in_coordinate[1]
                                 + first_normal_vector[2] * in_coordinate[2])
        self.surface_sign.append(second_normal_vector[0] * in_coordinate[0] + second_normal_vector[1] * in_coordinate[1]
                                 + second_normal_vector[2] * in_coordinate[2])
        surface[4].append(surface[4][0] * co_bottom[0] + surface[4][1] * co_bottom[1] + surface[4][2] * co_bottom[2])
        surface[5].append(surface[5][0] * co_bottom[0] + surface[5][1] * co_bottom[1] + surface[5][2] * co_bottom[2])
        for i in range(0, 6):
            self.surf_ids.append(max_surf_id + i + 1)
            surf = Surface.externalization(number=self.surf_ids[i], type=surf_type[i], parameters=surface[i],
                                           rotate=self.rotate, rotate_angle=self.rotate_angle, move=self.move)
            surf = surf.transfer()
            surface_model.append(surf)
        # 对圆柱面以及上下底面返回运算
        return Surfaces(surfaces=surface_model)

    def bound(self, in_=True):
        if not self.surf_ids:
            raise ValueError('surface id is empty!\n')
        if in_:
            bool_list = f'({self.surf_ids[0]}&{-self.surf_ids[1]}&{self.surf_ids[2]}&{-self.surf_ids[3]}'
            for i, value in enumerate(self.surface_sign):
                if value < 0:
                    bool_list += f'&{-self.surf_ids[4 + i]}'
                else:
                    bool_list += f'&{self.surf_ids[4 + i]}'
            bool_list += ')'
            return bool_list
        else:
            bool_list = f'({-self.surf_ids[0]}:{self.surf_ids[1]}:{-self.surf_ids[2]}:{self.surf_ids[3]}'
            for i, value in enumerate(self.surface_sign):
                if value < 0:
                    bool_list += f':{self.surf_ids[4 + i]}'
                else:
                    bool_list += f':{-self.surf_ids[4 + i]}'
            bool_list += ')'
            return bool_list


class Macrobodies(BaseModel):
    yaml_tag = u'!macrobodies'

    def __init__(self, macrobodies=None):
        self.macrobodies = macrobodies
        if self.macrobodies is None:
            self.macrobodies = []

    def check(self):
        pass

    def __str__(self):
        if len(self.macrobodies) > 0:
            s = 'MACROBODY\n'
            for body in self.macrobodies:
                s += str(body)
            s += '\n\n'
            return s
        else:
            return ''
