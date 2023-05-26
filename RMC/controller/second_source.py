# -*- coding:utf-8 -*-
# author: Hao Luo
# date: 2022-3-15
import h5py

from RMC.model.input.ExternalSource import *
from RMC.model.input.Physics import Physics
from RMC.parser.PlainParser import PlainParser


class SecondSource:
    """　由燃耗计算输出的二次源信息，生成固定源模式下的通用源
    Attributes
    __________
    step: int
        指定的燃耗步，在该燃耗步下，通过读取Result.h5文件，获得固定源模式下的通用源
    """
    def __init__(self, step=1):
        self._step = step

    @property
    def step(self):
        return self._step

    @step.setter
    def step(self, s):
        self._step = s

    def __call__(self, hdf5=None):
        """　解析燃耗计算生成的Result.h5文件，并生成固定源模式下可用的通用源分布(直接基于Cell抽样，注:仅支持简单几何)"""
        sources = []
        distributions = []
        result = h5py.File(hdf5, mode='r')

        if 'Step' + str(self._step) not in result['Burnup'].keys():
            raise AttributeError(f'Group Step{self._step} does not exist in hdf5 file {hdf5}\n')
        gamma_source = result['Burnup']['Step' + str(self._step)]['GammaSource']

        id = 0
        for cell in gamma_source.keys():
            id += 1

            sources.append(
                Source(source_id=id, fraction=1, particle=[2], energy=[f'd{id}'], cell=[cell]))

            energy = list(gamma_source[cell]['Energy'][...])
            intensity = list(gamma_source[cell]['Intensity'][...])

            distributions.append(Distribution(id=id, type=1, value=energy, probability=intensity))
        result.close()

        return ExternalSource(source=sources, distributions=distributions)
