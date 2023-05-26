from unittest import TestCase
from RMC.parser.PlainParser import PlainParser
import os

from pathlib import Path
dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()

def setUpModule():
    os.chdir(dir_path)

def tearDownModule():
    os.chdir(original_path)


class TestUniverse(TestCase):
    def test_compare_pattern(self):
        model_inp = 'resources/inp'

        geom = PlainParser(model_inp).parsed['geometry']
        core_univ = geom.get_univ(2)
        assem_univ1 = geom.get_univ(21)
        assem_univ2 = geom.get_univ(22)

        self.assertEqual([False, True], assem_univ1.compare_pattern(assem_univ2, excepts=[[1, 1]]))
        self.assertEqual([False, False], assem_univ1.compare_pattern(assem_univ2, excepts=[[1, 1], [2, 2]]))
        self.assertEqual([False, True], assem_univ2.compare_pattern(assem_univ1))
        self.assertRaises(ValueError, core_univ.compare_pattern, assem_univ2)
        self.assertRaises(ValueError, assem_univ2.compare_pattern, core_univ)
        self.assertRaisesRegex(ValueError, "different sizes", core_univ.compare_pattern, assem_univ2)

    def test_duplicate(self):
        model_inp = 'resources/inp'
        output_file = 'resources/output'
        reference = 'resources/reference_inp'

        model = PlainParser(model_inp).parsed
        geom = model['geometry']
        core_univ = geom.get_univ(2)
        assem1 = geom.get_univ(11)
        assem2 = geom.get_univ(12)
        assem_univ1 = geom.get_univ(21)
        assem_univ2 = geom.get_univ(22)

        new_univ1 = assem1.duplicate(geom, end_univ=assem_univ1, excepts=[[0, 0], [1, 1]],
                                     template_univ=assem_univ2)
        new_univ2 = assem1.duplicate(geom, end_univ=assem_univ1)
        new_univ3 = assem2.duplicate(geom, end_univ=assem_univ2, excepts=[[0, 0], [2, 2]],
                                     template_univ=assem_univ1)
        new_univ4 = assem2.duplicate(geom)

        core_univ.lattice.fill[3] = new_univ1[0].number
        core_univ.lattice.fill[5] = new_univ2[0].number
        core_univ.lattice.fill[1] = new_univ3[0].number
        core_univ.lattice.fill[4] = new_univ4[0].number

        with open(output_file, 'w') as f:
            f.write(str(model))

        output_model = PlainParser(output_file).parsed
        reference_model = PlainParser(reference).parsed

        self.assertEqual(str(output_model), str(reference_model))

        os.remove(output_file)
