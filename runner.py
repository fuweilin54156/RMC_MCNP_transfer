# -*- coding:utf-8 -*-
# author: Shen PF
# date: 2021-07-23
"""
examples:
>>> inp_MCNP = '06'
>>> M2R.transfer(inp_MCNP)

support * and ? in input (only one is supported):
'*' matches any char
'?' matches single char
for example:
'01?' matches ['011', '01a']
'01*' matches ['01', '011234567']
"""

import MCNPtoRMC as M2R
import os
import re

print("MCNP to RMC transformation tool V2.0 \nAuthor: Shen Pengfei\n")
# M2R.transfer('inp')
# files = ['0'+str(i+1) for i in range(8)]
# for file in files:
#     M2R.transfer(file)

filename = input(
    "Please input the filename: \nnote: '*' and '?' can be used once. " 
    "\n('*' can replace any char in any length; '?' can replace any single char)\n")
print("You have input : " + filename + '\n')
reg = filename.replace('?', '.')
reg = reg.replace('*', '.*')

processed_files = []
path = os.getcwd()
for file in os.listdir(path):
    if re.match(r'^' + reg + '$', file) and len(file) >= len(reg) - 1:
        processed_files.append(file)

print('These files are going to be processed: \n ' + ', '.join(processed_files) + '\n')
for file in processed_files:
    M2R.transfer(file)
