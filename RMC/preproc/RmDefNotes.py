import re
import os
import RMC.preproc.Base as FB


class RmDefNotes:
    def __init__(self, filename):
        self.file_name = os.path.abspath(filename)
        self.content = ''

    def rm_notes(self):
        self._read_in()
        '''
        删除变量定义语句思路：
        循环匹配，删除@开头的line
        '''
        note_find = re.search(r'[ ]*@.*\n', self.content)
        while note_find:
            indexs = note_find.span()
            self.content = FB.Base.replace_string(self.content, indexs, '')
            note_find = re.search(r'[ ]*@.*\n', self.content)

        self._write_file()

    def _read_in(self):
        with open(self.file_name, 'r') as f:
            self.content = f.read()

    def _write_file(self):
        with open(self.file_name, 'w+') as f:
            f.write(self.content)
