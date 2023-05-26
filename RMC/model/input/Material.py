# -*- coding: utf-8 -*-
# author: Xiaoyu Guo

from RMC.model.input.base import YMLModelObject as BaseModel


class Materials(BaseModel):
    card_option_types = {
        "MTLIB": {
            "PLIB": [str],
            "PNLIB": [str]
        },
        "CEACE": {
            "ERGBINHASH": [bool],
            "PTABLE": [bool],
            "OTFPTB": [bool],
            "OTFSAB": [bool],
            "DBRC": [bool],
            "TMS": [bool],
            "TMSTALLY": [bool],
            "OTFDB": [bool]
        },
        "MGACE": {
            "ERGGRP": [int]
        },
        "NUBAR": {
            "FACTOR": [int]
        },
    }

    def __init__(self, mats=None, mtlib=None, ceace=None, mgace=None, nubar=None, unparsed=''):
        if mats is None:
            mats = []
        self._mats = mats
        self._mtlib = mtlib
        self._ceace = ceace
        self._mgace = mgace
        self._nubar = nubar
        self._unparsed = unparsed  # ceace, mgace, etc.

    @property
    def mats(self):
        return self._mats

    def add_mat(self, mat):
        if mat not in self._mats:
            self._mats.append(mat)

    def update_mat(self, mat_id, nuclides):
        for mat in self._mats:
            if mat.mat_id == mat_id:
                for nuclide in nuclides:
                    mat.update_nuclide(nuclide)

    def check(self):
        pass

    def __str__(self):
        card = 'MATERIAL\n'
        for mat in self._mats:
            card += str(mat) + '\n'
        if self._mtlib:
            card += 'Mtlib'
            if "PLIB" in self._mtlib:
                card += ' plib=' + self._mtlib['PLIB']
            if "PNLIB" in self._mtlib:
                card += ' pnlib=' + self._mtlib['PNLIB']
            card += '\n'
        if self._ceace:
            card += "Ceace"
            if "ERGBINHASH" in self._ceace:
                card += ' ergbinhash=%d' % self._ceace['ERGBINHASH']
            if "PTABLE" in self._ceace:
                card += ' ptable=%d' % self._ceace['PTABLE']
            if "OTFPTB" in self._ceace:
                card += ' otfptb=%d' % self._ceace['OTFPTB']
            if "DBRC" in self._ceace:
                card += ' dbrc=%d' % self._ceace['DBRC']
            if "TMS" in self._ceace:
                card += ' tms=%d' % self._ceace['TMS']
            if "TMSTally" in self._ceace:
                card += ' tmstally=%d' % self._ceace['TMSTALLY']
            if "OTFDB" in self._ceace:
                card += ' otfdb=%d' % self._ceace['OTFDB']
            card += '\n'
        if self._mgace:
            card += "Mgace "
            card += "erggrp=" + str(self._mgace['ERGGRP'])
            card += '\n'
        if self._nubar:
            card += "Nubar "
            card += 'factor=' + str(self._nubar["FACTOR"])
            card += '\n'
        if self._unparsed != '':
            card += self._unparsed
        card += '\n\n'
        return card


class Material(BaseModel):
    def __init__(self, mat_id=0, density=0, nuclides=None):
        self._mat_id = mat_id
        self._density = density

        if nuclides is None:
            nuclides = []
        self._nuclides = nuclides

    @property
    def mat_id(self):
        return self._mat_id

    @mat_id.setter
    def mat_id(self, mat_id):
        self._mat_id = mat_id

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, density):
        self._density = density

    def check(self):
        assert isinstance(self._mat_id, int)

    def update_nuclide(self, nuclide):
        for current_nuclide in self._nuclides:
            if current_nuclide.name == nuclide.name:
                current_nuclide.density = nuclide.density
                return

        self._nuclides.append(nuclide)

    def __str__(self):
        card = 'Mat'
        card += ' ' + str(self._mat_id)
        card += ' ' + str(self._density)
        for nuclide in self._nuclides:
            card += '\n ' + str(nuclide)
        return card


class SabMaterial(BaseModel):
    def __init__(self, mat_id=0, sab_nuclides=None):
        self._mat_id = mat_id
        self._nuclides = sab_nuclides

    @property
    def mat_id(self):
        return self._mat_id

    @property
    def nuclides(self):
        return self._nuclides

    def __str__(self):
        card = 'Sab '
        card += str(self._mat_id) + ' '
        card += ' '.join(self._nuclides)
        return card


class DynamicMaterial(BaseModel):
    card_option_types = {
        "TIME": ['list', float, -1],
        "MATDENVALUE": ['list', float, -1],
        "NUCDENVALUE": ['list', float, -1]
    }

    def __init__(self, mat_id=0, options=None):
        self._mat_id = mat_id
        self._options = options

    @property
    def mat_id(self):
        return self._mat_id

    @property
    def options(self):
        return self._options

    def __str__(self):
        card = 'DynamicMat %d' % self._mat_id
        if "TIME" in self._options:
            card += " time = " + ' '.join([str(x) for x in self._options['TIME']])
        if "MATDENVALUE" in self._options:
            card += " matdenvalue = " + ' '.join([str(x) for x in self._options['MATDENVALUE']])
        if "NUCDENVALUE" in self._options:
            card += " nucdenvalue = " + ' '.join([str(x) for x in self._options['NUCDENVALUE']])
        return card


class Nuclide(BaseModel):
    def __init__(self, name=None, density=0):
        self._name = name
        self._density = density

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, density):
        self._density = density

    def check(self):
        pass

    def __str__(self):
        return self._name + ' ' + str(self._density)
