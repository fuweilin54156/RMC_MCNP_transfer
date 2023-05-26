#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2020-06-10


def dict2list(dic, null=None):
    result = []
    for key in dic:
        if dic[key] is None:
            result.append([key, null])
        else:
            result.append([key, dic[key]])
    return result
