# -*- codint:utf-8 -*-

import os
import shutil
import unittest
from unittest import TestCase
from RMC.parser.PlainParser import PlainParser

from pathlib import Path
dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()

def setUpModule():
    os.chdir(dir_path)

def tearDownModule():
    os.chdir(original_path)


class TestModel(TestCase):
    model = None
    base_dir = None
    all_original_file = None

    @classmethod
    def setUp(cls):
        cls.base_dir = os.getcwd()
        cls.all_original_file = os.listdir(cls.base_dir)

    @classmethod
    def tearDown(cls):
        all_file = os.listdir(cls.base_dir)
        all_new_file = list(set(all_file) - set(cls.all_original_file))
        for file in all_new_file:
            os.remove(os.path.join(cls.base_dir, file))

    @staticmethod
    def read_in(filename):
        with open(filename, 'r+') as f:
            content = f.read()
        return content

    def test_lat5(self):
        inp = 'resources/test_lat5/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_lat5/reference')
        self.assertEqual(str(self.model), test_result)

    def test_body(self):
        inp = 'resources/test_body/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_body/reference')
        self.assertEqual(str(self.model), test_result)

    def test_burnup(self):
        inp = 'resources/test_burnup/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_burnup/reference')
        self.assertEqual(str(self.model), test_result)

    def test_external_source(self):
        inp = 'resources/test_external_source/inp'
        shutil.copyfile(inp, 'inp')
        shutil.copyfile('resources/test_external_source/result.h5', 'result.h5')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_external_source/reference')
        self.assertEqual(str(self.model), test_result)

    def test_fixed_source(self):
        inp = 'resources/test_fixed_source/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_fixed_source/reference')
        self.assertEqual(str(self.model), test_result)

    def test_geometry(self):
        inp = 'resources/test_geometry/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_geometry/reference')
        self.assertEqual(str(self.model), test_result)

    def test_tally(self):
        inp = 'resources/test_tally/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_tally/reference')
        self.assertEqual(str(self.model), test_result)

    def test_material(self):
        inp = 'resources/test_material/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_material/reference')
        self.assertEqual(str(self.model), test_result)

    def test_surface_transformation(self):
        inp = 'resources/test_surface_transformation/inp'
        shutil.copyfile(inp, 'inp')
        self.model = PlainParser('inp').parsed
        test_result = TestModel.read_in('resources/test_surface_transformation/reference')
        self.assertEqual(str(self.model), test_result)


if __name__ == '__main__':
    unittest.main()
