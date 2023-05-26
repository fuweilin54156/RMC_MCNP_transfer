# -*- coding:utf-8 -*-
# author: Hao Luo
# date: 2022-3-15

from RMC.model.input.base import YMLModelObject as BaseModel


class Physics(BaseModel):
    card_option_types = {
        'PARTICLEMODE': ['list', str, -1],
        'NEUTRON': {
            'MINENERGY': [float]
        },
        'PHOTON': {
            'PPRODUCEE': [bool],
            'TTB': [bool],
            'ANNIHILATION': [bool],
            'COHERENT': [bool],
            'PHOTONNUCLEUS': [bool],
            'DOPPLER': [bool],
            'UPERERG': [float],
            'ERGCUTGMA': [float],
            'ELECMULTITIMES': [float]
        },
        'ELECTRON': {
            'MAXENERGY': [float],
            'MINENERGY': [float],
            'EPRODUCEP': [bool],
            'ERGLOSSSTRAGGLE': [bool],
            'BREMS': [bool],
            'BREMSANGLE': [bool],
            'BREMSPHOTONMULTITIMES': [float],
            'BREMSEACHSUBTEP': [bool],
            'BREMSERGLOSSMETHOD': [int],
            'XRAYMULTITIMES': [float],
            'KNOCKONMULTITIMES': [float]
        }
    }

    def __init__(self, particle_mode=None, neutron=None, photon=None, electron=None):
        self._particle_mode = particle_mode
        self._neutron = neutron
        self._photon = photon
        self._electron = electron

    @property
    def particle_mode(self):
        return self._particle_mode

    @property
    def neutron(self):
        return self._neutron

    @property
    def photon(self):
        return self._photon

    @property
    def electron(self):
        return self._electron

    def check(self):
        if self._electron and 'BREMSEACHSUBTEP' in self._electron:
            if self._electron['BREMSEACHSUBTEP'] and self._electron['BREMSPHOTONMULTITIMES'] != 1:
                raise ValueError("wrong BREMSPHOTONMULTITIMES value when BREMSEACHSUBTEP == 1")

    def __str__(self):
        card = 'PHYSICS\n'
        card += 'ParticleMode ' + ' '.join(self._particle_mode) + '\n'
        if self._neutron is not None:
            card += 'Neutron minenergy = ' + str(self._neutron['MINENERGY'])
            card += '\n'
        if self._photon is not None:
            card += 'Photon'
            if 'PPRODUCEE' in self._photon:
                card += ' pproducee = ' + '{:d}'.format(self._photon['PPRODUCEE'])
            if 'TTB' in self._photon:
                card += ' ttb = ' + '{:d}'.format(self._photon['TTB'])
            if 'ANNIHILATION' in self._photon:
                card += ' annihilation = ' + '{:d}'.format(self._photon['ANNIHILATION'])
            if 'COHERENT' in self._photon:
                card += ' coherent = ' + '{:d}'.format(self._photon['COHERENT'])
            if 'PHOTONNUCLEUS' in self._photon:
                card += ' photonnucleus = ' + '{:d}'.format(self._photon['PHOTONNUCLEUS'])
            if 'DOPPLER' in self._photon:
                card += ' doppler = ' + '{:d}'.format(self._photon['DOPPLER'])
            if 'UPERERG' in self._photon:
                card += ' upererg = ' + str(self._photon['UPERERG'])
            if 'ERGCUTGMA' in self._photon:
                card += ' ergcutgma = ' + str(self._photon['ERGCUTGMA'])
            if 'ELECMULTITIMES' in self._photon:
                card += ' elecmultitimes = ' + str(self._photon['ELECMULTITIMES'])
            card += '\n'
        if self._electron is not None:
            card += 'Electron'
            if 'MAXENERGY' in self._electron:
                card += ' maxenergy = ' + str(self._electron['MAXENERGY'])
            if 'MINENERGY' in self._electron:
                card += ' minenergy = ' + str(self._electron['MINENERGY'])
            if 'EPRODUCEP' in self._electron:
                card += ' eproducep = ' + '{:d}'.format(self._electron['EPRODUCEP'])
            if 'ERGLOSSSTRAGGLE' in self._electron:
                card += ' erglossstraggle = ' + '{:d}'.format(self._electron['ERGLOSSSTRAGGLE'])
            if 'BREMS' in self._electron:
                card += ' brems = ' + '{:d}'.format(self._electron['BREMS'])
            if 'BREMSANGLE' in self._electron:
                card += ' bremsangle = ' + '{:d}'.format(self._electron['BREMSANGLE'])
            if 'BREMSPHOTONMULTITIMES' in self._electron:
                card += ' bremsphotonmultitimes = ' + str(self._electron['BREMSPHOTONMULTITIMES'])
            if 'BREMSEACHSUBTEP' in self._electron:
                card += ' bremseachsubtep = ' + '{:d}'.format(self._electron['BREMSEACHSUBTEP'])
            if 'BREMSERGLOSSMETHOD' in self._electron:
                card += ' bremserglossmethod = ' + str(self._electron['BREMSERGLOSSMETHOD'])
            if 'XRAYMULTITIMES' in self._electron:
                card += ' xraymultitimes = ' + str(self._electron['XRAYMULTITIMES'])
            if 'KNOCKONMULTITIMES' in self._electron:
                card += ' knockonmultitimes = ' + str(self._electron['KNOCKONMULTITIMES'])
            card += '\n'
        card += '\n\n'
        return card
