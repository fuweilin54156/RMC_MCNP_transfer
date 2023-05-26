from RMC.model.input.base import YMLModelObject as BaseModel
from RMC.model.input.Surface import *


class Binaryout(BaseModel):
    yaml_tag = u'!binaryout'
    card_option_types = {
        'WRITE': [int],
        'SYMM': [int],
        'PARTYPE': [int],
        'SURF': ['list', float],
        'CEL': ['list', int],
        'WCYCLES': [int]
    }

    def __init__(self, write=None, symm=None, partype=None, surf=None, cel=None, wcycles=None):
        self._write = write
        self._symm = symm
        self._partype = partype
        self._surf = surf
        self._cel = cel
        self._wcycles = wcycles

    def check(self):
        pass

    @property
    def surf(self):
        return self._surf

    @surf.setter
    def surf(self, surf):
        self._surf = surf

    def pro_binaryout(self, macrobodies=None):
        self.surf = self.surface(self.surf, macrobodies)

    def __str__(self):
        card = 'BINARYOUT\nWrtSurfSrc '
        if self._write is not None:
            card += ' Write = ' + str(self._write)
        if self._symm is not None:
            card += ' SYMM = ' + str(self._symm)
        if self._partype is not None:
            card += ' PARTYPE = ' + str(self._partype)
        if self._surf is not None:
            card += ' Surf = ' + ' '.join([str(int(x)) if isinstance(x, float) else str(x) for x in self._surf])
        if self._cel is not None:
            card += ' Cel = ' + ' '.join([str(x) for x in self._cel])
        if self._wcycles is not None:
            card += ' WCycles = ' + str(self._wcycles)
        card += '\n\n'
        return card

    def add_options(self, options):
        if 'WRITE' in options.keys():
            self._write = options['WRITE']
        if 'SYMM' in options.keys():
            self._symm = options['SYMM']
        if 'PARTYPE' in options.keys():
            self._partype = options['PARTYPE']
        if 'SURF' in options.keys():
            self._surf = options['SURF']
        if 'CEL' in options.keys():
            self._cel = options['CEL']
        if 'WCYCLES' in options.keys():
            self._wcycles = options['WCYCLES']

    @staticmethod
    def surface(surf_list, macrobodies):
        flag_str = False
        for i, surface in enumerate(surf_list):
            if surface == '(':
                flag_str = True
                continue
            elif surface == ')':
                flag_str = False
                continue
            elif flag_str:
                continue
            else:
                if abs(surface) - int(abs(surface)) == 0:
                    if int(abs(surface)) in macrobodies.keys():
                        surf_list[i] = macrobodies[int(abs(surface))].surf_ids[0]
                    else:
                        continue
                else:
                    surf_list[i] = macrobodies[int(abs(surface))].surf_ids[int(str(surface).split('.')[1]) - 1]
                if surface < 0:
                    surf_list[i] = -surf_list[i]
        return surf_list
