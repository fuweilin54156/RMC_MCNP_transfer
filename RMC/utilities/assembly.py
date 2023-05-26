#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2020-06-10

import re


class AssemblyRecord:
    curCycleType = "cur"
    skipCycleType = "skip"
    newAssemType = "new"

    def __init__(self, pos: int, meta: str, assem: int):
        self._pos = pos
        self._meta = meta
        self._assem = assem

    @property
    def pos(self):
        return self._pos

    @property
    def meta(self):
        return self._meta

    @property
    def assem(self):
        return self._assem

    @pos.setter
    def pos(self, pos):
        self._pos = pos

    @meta.setter
    def meta(self, meta):
        self._meta = meta

    @assem.setter
    def assem(self, assem):
        self._assem = assem


def extract_alias(alias, univ_coordinate):
    pattern = re.compile(r'([A-Za-z_]+)(\d+)(?:-(\d+))?')
    m = pattern.match(univ_coordinate)
    n_column = len(alias['column'])
    if m:
        if m.group(1) == alias['new']:
            assem = int(m.group(2))
            return AssemblyRecord(pos=-1, meta=AssemblyRecord.newAssemType, assem=assem)
        else:
            pos = n_column * alias['row'].index(int(m.group(2))) \
                  + alias['column'].index(m.group(1))
            if m.group(3) is None:
                return AssemblyRecord(pos=pos, meta=AssemblyRecord.curCycleType, assem=0)
            else:
                cycle = int(m.group(3))
                return AssemblyRecord(pos=pos, meta=AssemblyRecord.skipCycleType, assem=-cycle)

    else:
        raise ValueError('Entry %s in refuelling can not be recognized!' % univ_coordinate)
