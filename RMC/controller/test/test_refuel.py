from unittest import TestCase
from RMC.controller.refuel import Refuel
from RMC.parser.PlainParser import PlainParser
import filecmp
import os
import shutil
import numpy as np

import unittest

from pathlib import Path
dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()

def setUpModule():
    os.chdir(dir_path)

def tearDownModule():
    os.chdir(original_path)


class TestRefuel(TestCase):
    def test_refuel_control_rod2(self):
        """

        :return:
        """
        # input files
        inp = 'resources/2/inp'
        inp_initial = 'resources/2/inp_initial'
        refuel_inp = 'resources/2/refuelling.yml'
        status_file = 'resources/2/status.txt'
        mat_ref = ['resources/2/mat_1_ref.npy', 'resources/2/mat_2_ref.npy']
        mat_bak = ['resources/2/mat_1_bak.npy', 'resources/2/mat_2_bak.npy']
        # material file will be changed in the test, so the mat_1.npy file is copied from mat_1_bak.npy
        #   in each time of testing.
        mat = ['resources/2/mat_1.npy', 'resources/2/mat_2.npy']

        # output files
        output = 'resources/2/inp.refuel'

        # reference
        reference = 'resources/2/reference'

        # storage
        storage = 'resources/2/storage'

        for idx in range(len(mat_bak)):
            shutil.copy(mat_bak[idx], mat[idx])

        # 0. Load the data and do a step of refuelling.
        refuel = Refuel(refuel_inp=refuel_inp)
        base_model = PlainParser(inp_initial).parsed
        refuel.refuel(1, inp, base_model, output=output, save_remove={
            "storage": storage,
            "cycle": 1,
        })

        # 1. Check the plain input file
        self.assertTrue(filecmp.cmp(output, reference))
        # assert str(PlainParser(output).parsed) == str(PlainParser(reference).parsed)

        # 2. Check the material npy file
        for idx in range(len(mat)):
            a = np.load(mat[idx])
            b = np.load(mat_ref[idx])
            try:
                self.assertTrue(np.all(a == b))
            except AssertionError:
                print("output {}: ".format(idx + 1))
                print(a)
                print("reference {}: ".format(idx + 1))
                print(b)
                exit(1)

        # 3. Postprocessing
        os.remove(output)
        for idx in range(len(mat)):
            os.remove(mat[idx])

        shutil.rmtree(storage)

    def test_refuel_control_rod3(self):
        """

        :return:
        """
        # input files
        folders = ['3', '4', '5', '6', '7', '8', '9']
        mat_nums = [2, 2, 4, 4, 4, 4, 4]
        steps = [3, 3, 1, 3, 1, 3, 3]
        cycles = [2, 2, 1, 2, 1, 2, 2]
        # folders = ['8']
        # mat_nums = [4]
        # steps = [3]
        for idx in range(len(folders)):
            folder = folders[idx]
            inp = 'resources/{}/inp'.format(folder)
            inp_initial = 'resources/{}/inp_initial'.format(folder)
            refuel_inp = 'resources/{}/refuelling.yml'.format(folder)
            refuel_inp_bak = 'resources/{}/refuelling_bak.yml'.format(folder)
            refuel_inp_ref = 'resources/{}/refuelling_ref.yml'.format(folder)
            status_file = 'resources/{}/status.txt'.format(folder)
            mat_ref = ['resources/{}/mat_{}_ref.npy'.format(folder, str(i)) for i in range(1, mat_nums[idx] + 1)]
            mat_bak = ['resources/{}/mat_{}_bak.npy'.format(folder, str(i)) for i in range(1, mat_nums[idx] + 1)]
            # material file will be changed in the test, so the mat_1.npy file is copied from mat_1_bak.npy
            #   in each time of testing.
            mat = ['resources/{}/mat_{}.npy'.format(folder, str(i)) for i in range(1, mat_nums[idx] + 1)]

            # output files
            output = 'resources/{}/inp.refuel'.format(folder)

            # reference
            reference = 'resources/{}/reference'.format(folder)

            # storage
            storage = 'resources/{}/storage'.format(folder)
            storage_bak = 'resources/{}/storage_bak'.format(folder)
            # todo: check the output of assembly storage.
            #       currently only assemblies from cycle 1 can be fetched.
            if os.path.exists(os.path.join(storage, "cycle1")):
                shutil.rmtree(os.path.join(storage, "cycle1"))
            if os.path.exists(os.path.join(storage_bak, "cycle1")):
                shutil.copytree(os.path.join(storage_bak, "cycle1"), os.path.join(storage, "cycle1"))

            for mat_idx in range(len(mat_bak)):
                shutil.copy(mat_bak[mat_idx], mat[mat_idx])
            shutil.copy(refuel_inp_bak, refuel_inp)

            # 0. Load the data and do a step of refuelling.
            refuel = Refuel(refuel_inp=refuel_inp)
            base_model = PlainParser(inp_initial).parsed
            refuel.refuel(steps[idx], inp, base_model, output=output, dump=True, save_remove={
                "storage": storage,
                "cycle": cycles[idx],
            })

            # 1. Check the plain input file
            a = PlainParser(reference).parsed
            b = PlainParser(output).parsed

            self.assertEqual(str(a), str(b))
            # self.assertTrue(filecmp.cmp(output, reference))
            self.assertTrue(filecmp.cmp(refuel_inp, refuel_inp_ref))

            # 2. Check the material npy file
            for mat_idx in range(len(mat)):
                a = np.load(mat[mat_idx])
                b = np.load(mat_ref[mat_idx])
                try:
                    self.assertTrue(np.all(a == b))
                except AssertionError:
                    print("output {}: ".format(mat_idx + 1))
                    print(a)
                    print("reference {}: ".format(mat_idx + 1))
                    print(b)
                    exit(1)

            # 3. Postprocessing
            os.remove(output)
            for mat_idx in range(len(mat)):
                os.remove(mat[mat_idx])
            os.remove(refuel_inp)
            shutil.rmtree(storage)
            print("Test in Folder {} passed".format(folder))

    # todo: add a test case with unused universe in the lattice.

    def test_refuel_status(self):
        """

        :return:
        """
        # input files
        inp_initial = 'resources/1/inp_initial'
        refuel_inp = 'resources/1/refuelling.yml'
        inp = 'resources/1/inp'
        status_file = 'resources/1/status.txt'
        mat_ref = 'resources/1/mat_1_ref.npy'
        mat_bak = 'resources/1/mat_1_bak.npy'
        # material file will be changed in the test, so the mat_1.npy file is copied from mat_1_bak.npy
        #   in each time of testing.
        mat = 'resources/1/mat_1.npy'

        # output files
        output = 'resources/1/inp.refuel'

        # reference
        reference = 'resources/1/reference'

        # storage
        storage = 'resources/1/storage'

        shutil.copy(mat_bak, mat)

        # 0. Load the data and do a step of refuelling.
        refuel = Refuel(refuel_inp=refuel_inp)
        base_model = PlainParser(inp_initial).parsed
        refuel.refuel(1, inp, base_model, output=output, save_remove={
            "storage": storage,
            "cycle": 1,
        })

        # 1. Check the plain input file
        self.assertTrue(filecmp.cmp(output, reference))
        # assert str(PlainParser(output).parsed) == str(PlainParser(reference).parsed)

        # 2. Check the material npy file
        a = np.load(mat)
        b = np.load(mat_ref)
        try:
            self.assertTrue(np.all(a == b))
        except AssertionError:
            print("output: ")
            print(a)
            print("reference: ")
            print(b)
            exit(1)

        # 3. Postprocessing
        os.remove(output)
        os.remove(mat)
        shutil.rmtree(storage)

    # todo: add a test case with unused universe in the lattice.


if __name__ == '__main__':
    unittest.main()
