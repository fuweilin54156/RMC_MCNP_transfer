# -*- codint:utf-8 -*-

import os

from RMC.controller.ctf import FakeCTF
from RMC.controller.ctf import FakeCTFPreproc

import unittest
from unittest import TestCase


from pathlib import Path
dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()

def setUpModule():
    os.chdir(dir_path)

def tearDownModule():
    os.chdir(original_path)


class TestFakeCTF(TestCase):
    fake_ctf = None
    fake_ctf_preproc = None
    base_dir = None

    @classmethod
    def setUp(cls):
        inp = 'resources/FakeCTF/inp'
        cls.fake_ctf = FakeCTF(inp)
        cls.fake_ctf_preproc = FakeCTFPreproc(inp)
        cls.base_dir = os.path.dirname(inp)

    @classmethod
    def tearDown(cls):
        os.remove(os.path.join(cls.base_dir, 'deck.ctf.h5'))

    def test_run(self):
        self.fake_ctf_preproc.run()
        self.fake_ctf.run()
        self.assertTrue(os.path.exists(os.path.join(self.base_dir, 'deck.ctf.h5')))


class TestFakeCTFPreproc(TestCase):
    fake_ctf_preproc = None
    base_dir = None
    original_path = None

    @classmethod
    def setUp(cls):
        inp = 'resources/FakeCTF/inp'
        cls.fake_ctf_preproc = FakeCTFPreproc(inp)
        cls.base_dir = os.path.dirname(inp)

        TestFakeCTF.original_path = os.getcwd()
        os.chdir(dir_path)

    @classmethod
    def tearDown(cls):
        os.remove(os.path.join(cls.base_dir, 'deck.inp'))

        os.chdir(TestFakeCTF.original_path)

    def test_run(self):
        self.fake_ctf_preproc.run()
        self.assertTrue(os.path.exists(os.path.join(self.base_dir, 'deck.inp')))


if __name__ == '__main__':
    unittest.main()