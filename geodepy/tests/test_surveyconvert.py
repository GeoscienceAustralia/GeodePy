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
        files = ['fbksample.fbk', 'fbksample.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_fbk_filepath = os.path.join(tempdir.name, 'fbksample.fbk')
        temp_fbkgpy_filepath = os.path.join(tempdir.name, 'fbksample.gpy')
        fbk2msr(temp_fbk_filepath, temp_fbkgpy_filepath, strict=False, zerodist=False, same_stdev=False)
        original_msr = open(os.path.join(abs_path, 'resources/fbksample.msr'))
        temp_msr = open(os.path.join(tempdir.name, 'fbksample.msr'))
        assert [row for row in original_msr] == [row for row in temp_msr]
        original_msr.close()
        temp_msr.close()
        tempdir.cleanup()

    def test_writestn(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['fbksample.txt', 'fbksample.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_txt_filepath = os.path.join(tempdir.name, 'fbksample.txt')
        writestn(temp_txt_filepath, 'S56')
        original_stn = open(os.path.join(abs_path, 'resources/fbksample.stn'))
        temp_stn = open(os.path.join(tempdir.name, 'fbksample.stn'))
        assert [row for row in original_stn] == [row for row in temp_stn]
        original_stn.close()
        temp_stn.close()
        tempdir.cleanup()

    def test_gsi2msr(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['gsisample.gsi', 'gsisample.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_gsi_filepath = os.path.join(tempdir.name, 'gsisample.gsi')
        temp_gsigpy_filepath = os.path.join(tempdir.name, 'gsisample.gpy')
        gsi2msr(temp_gsi_filepath, temp_gsigpy_filepath)
        original_msr = open(os.path.join(abs_path, 'resources/gsisample.msr'))
        temp_msr = open(os.path.join(tempdir.name, 'gsisample.msr'))
        assert [row for row in original_msr] == [row for row in temp_msr]
        original_msr.close()
        temp_msr.close()
        tempdir.cleanup()

    def test_gsi2msr_no_config(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['gsisample.gsi']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_gsi_filepath = os.path.join(tempdir.name, 'gsisample.gsi')
        gsi2msr(temp_gsi_filepath)
        original_stn = open(os.path.join(abs_path, 'resources/gsisample_noconfig.msr'))
        temp_stn = open(os.path.join(tempdir.name, 'gsisample.msr'))
        assert [row for row in original_stn] == [row for row in temp_stn]
        original_stn.close()
        temp_stn.close()
        tempdir.cleanup()

    def test_gsi2stn(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['gsisample.gsi']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_gsi_filepath = os.path.join(tempdir.name, 'gsisample.gsi')
        #TODO: finish test_gsi2stn testing
        tempdir.cleanup()
        pass


if __name__ == '__main__':
    unittest.main()
