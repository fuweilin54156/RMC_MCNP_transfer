# -*- codint:utf-8 -*-

import os
import filecmp
import shutil

import unittest
from unittest import TestCase
from RMC.controller.rmc import RMCController, FakeRMC
from RMC.controller.second_source import SecondSource

from pathlib import Path
dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()

def setUpModule():
    os.chdir(dir_path)

def tearDownModule():
    os.chdir(original_path)


class TestSecondSource(TestCase):
    controller = None
    base_dir = None
    all_original_file = None

    @classmethod
    def setUp(cls):
        inp = 'resources/FakeSecondSource/inp'
        cls.base_dir = os.path.dirname(inp)
        archive = os.path.join(cls.base_dir, 'archive')
        cls.controller = RMCController(inp, archive)
        cls.all_original_file = os.listdir(cls.base_dir)

    @classmethod
    def tearDown(cls):
        all_file = os.listdir(cls.base_dir)
        all_new_file = list(set(all_file) - set(cls.all_original_file))
        for file_dir in all_new_file:
            file_dir = os.path.join(cls.base_dir, file_dir)
            if os.path.isfile(file_dir):
                os.remove(file_dir)
            else:
                shutil.rmtree(file_dir)

    def test_second_source(self):
        all_original_file = os.listdir(self.base_dir)
        controller_property = {'inp': 'inp', 'archive': 'archive'}
        while self.controller.continuing('resources/status.txt', controller_property):
            fake_rmc = FakeRMC(controller_property['inp'], controller_property['archive_dir'])
            fake_rmc.run()
            fake_rmc.archive_output()
            shutil.copy("resources/FakeSecondSource/result.h5",
                        "resources/FakeSecondSource/archive/burnup/inp.burnup.Result.h5")

        all_file = os.listdir(self.base_dir)
        all_new_file = list(set(all_file) - set(all_original_file))
        for file in all_new_file:
            if not os.path.isfile(file):
                continue
            self.assertTrue(
                filecmp.cmp(os.path.join(self.base_dir, file),
                            os.path.join(self.base_dir, 'reference', file)),
                'file {} is not matched!'.format(file))