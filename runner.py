# -*- coding:utf-8 -*-
# author: Shen PF
# date: 2021-07-23
"""
examples:
>>> import MCNPtoRMC as M2R
>>> inp_MCNP = '06'
>>> M2R.transfer(inp_MCNP)

support * and ? in input (only one is supported):
'*' matches any char
'?' matche
s single char
for example:
'01?' matches ['011', '01a']
'01*' matches ['01', '011234567']
"""

import MCNPtoRMC as M2R
import os
import re

print('-------------------------------------------------------------------------\n'
      "                  MM   MM      2222      RRRRR \n"
      "                  M M M M          2     RR   R \n"
      "                  M  M  M       22       RR RR  \n"
      "                  M     R     222222     RR   R \n\n"
      "MCNP to RMC transformation tool V2.0 \nAuthor: Shen Pengfei, Gou Yuanhao\n"
      "Email: 2043965149@qq.com\n"
      '-------------------------------------------------------------------------\n')
# M2R.transfer('fq-HZP-MCNP')
# files = ['0'+str(i+1) for i in range(8)]
# for file in files:
#     M2R.transfer(file)

filename = input(
    "Please input the filename: \nnote: '*' and '?' can be used once. " 
    "'*' can replace any char in any length; '?' can replace any single char.\n"
    "for example: '01?' matches ['011', '01a'] ;\n"
    "             '01*' matches ['01', '011234567']\n"
    ">> ")
# print("You have input : " + filename + '\n')
reg = filename.replace('?', '.')
reg = reg.replace('*', '.*')

processed_files = []
path = os.getcwd()
for file in os.listdir(path):
    if re.match(r'^' + reg + '$', file) and len(file) >= len(reg) - 1:
        processed_files.append(file)

print('\n-------------------------------------------------------------------------\n'
      'These files are going to be processed: \n ' + ', '.join(processed_files) + '\n'
      '-------------------------------------------------------------------------\n')
for file in processed_files:
    M2R.transfer(file)

os.system("pause")
