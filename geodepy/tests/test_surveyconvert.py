import unittest
import tempfile
import shutil
import os.path
from geodepy.surveyconvert.fbk import *
from geodepy.surveyconvert.gsi import *


class TestSurveyConvert(unittest.TestCase):
    def test_fbk2msr(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['Site13-134.fbk', 'Site13-134.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_fbk_filepath = os.path.join(tempdir.name, 'Site13-134.fbk')
        temp_fbkgpy_filepath = os.path.join(tempdir.name, 'Site13-134.gpy')
        fbk2msr(temp_fbk_filepath, temp_fbkgpy_filepath, strict=False, zerodist=False, same_stdev=False)
        original_msr = open(os.path.join(abs_path, 'resources/Site13-134.msr'))
        temp_msr = open(os.path.join(tempdir.name, 'Site13-134.msr'))
        assert [row for row in original_msr] == [row for row in temp_msr]
        original_msr.close()
        temp_msr.close()

    def test_writestn(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['Site13-134.txt', 'Site13-134.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_txt_filepath = os.path.join(tempdir.name, 'Site13-134.txt')
        writestn(temp_txt_filepath, 'S56')
        original_stn = open(os.path.join(abs_path, 'resources/Site13-134.stn'))
        temp_stn = open(os.path.join(tempdir.name, 'Site13-134.stn'))
        assert [row for row in original_stn] == [row for row in temp_stn]
        original_stn.close()
        temp_stn.close()

    def test_gsi2msr(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['ST0618HZ.gsi', 'ST0618HZ.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_gsi_filepath = os.path.join(tempdir.name, 'ST0618HZ.gsi')
        temp_gsigpy_filepath = os.path.join(tempdir.name, 'ST0618HZ.gpy')
        gsi2msr(temp_gsi_filepath, temp_gsigpy_filepath)
        original_stn = open(os.path.join(abs_path, 'resources/ST0618HZ.msr'))
        temp_stn = open(os.path.join(tempdir.name, 'ST0618HZ.msr'))
        assert [row for row in original_stn] == [row for row in temp_stn]
        original_stn.close()
        temp_stn.close()


if __name__ == '__main__':
    unittest.main()