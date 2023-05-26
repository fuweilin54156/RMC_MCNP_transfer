import unittest
import filecmp
import os
from RMC.modifier.formatter.PlainFormatter import PlainFormatter

from pathlib import Path
dir_path = Path(os.path.dirname(os.path.abspath(__file__)))

original_path = os.getcwd()

def setUpModule():
    os.chdir(dir_path)

def tearDownModule():
    os.chdir(original_path)


class TestPlainFormatter(unittest.TestCase):
    file_name = ['resources/BEAVRS/plain_inp', 'resources/assembly/plain_inp']

    def test_plain_formatter(self):
        output_file_name = ['resources/BEAVRS/formatted_plain_inp', 'resources/assembly/formatted_plain_inp']
        reference = ['resources/BEAVRS/reference_plain_inp', 'resources/assembly/reference_plain_inp']

        for i in range(len(self.file_name)):
            formatter = PlainFormatter(self.file_name[i])
            formatter.format_and_output(output_file_name[i])
            self.assertTrue(filecmp.cmp(output_file_name[i], reference[i]))
            os.remove(output_file_name[i])
            formatter.clear()

    def test_plain_format_to_card_list(self):
        output_file_name = ['resources/BEAVRS/formatted_card_list_plain_inp',
                            'resources/assembly/formatted_card_list_plain_inp']
        reference = ['resources/BEAVRS/reference_card_list_plain_inp',
                     'resources/assembly/reference_card_list_plain_inp']

        for i in range(len(self.file_name)):
            formatter = PlainFormatter(self.file_name[i])
            data = formatter.format_to_cards()
            with open(output_file_name[i], 'w') as f:
                f.write('\n\n'.join(data))

            self.assertTrue(filecmp.cmp(output_file_name[i], reference[i]))
            os.remove(output_file_name[i])
            formatter.clear()


if __name__ == '__main__':
    unittest.main()
