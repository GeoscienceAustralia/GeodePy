import unittest
import tempfile
import shutil
import os.path
from geodepy.surveyconvert.fbk import *


class TestSurveyConvert(unittest.TestCase):
    def test_fbk2msr(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['Site13-134.fbk', 'Site13-134.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_fbk_filepath = os.path.join(tempdir.name, 'Site13-134.fbk')
        temp_gpy_filepath = os.path.join(tempdir.name, 'Site13-134.gpy')
        fbk2msr(temp_fbk_filepath, temp_gpy_filepath, strict=False, zerodist=False, same_stdev=False)
        assert ([row for row in open(os.path.join(abs_path, 'resources/Site13-134.msr'))]
                == [row for row in open(os.path.join(tempdir.name, 'Site13-134.msr'))])

    def test_writestn(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        tempdir = tempfile.TemporaryDirectory()
        files = ['Site13-134.txt', 'Site13-134.gpy']
        for f in files:
            shutil.copy(os.path.join(abs_path, 'resources', f), tempdir.name)
        temp_txt_filepath = os.path.join(tempdir.name, 'Site13-134.txt')
        writestn(temp_txt_filepath, 'S56')
        assert([row for row in open(os.path.join(abs_path, 'resources/Site13-134.stn'))]
               == [row for row in open(os.path.join(tempdir.name, 'Site13-134.stn'))])


if __name__ == '__main__':
    unittest.main()