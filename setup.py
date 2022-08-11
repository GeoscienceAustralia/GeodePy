from setuptools import setup

setup(name='geodepy',
      version='0.4.0',
      description='GA Geodesy Package',
      long_description='A toolkit for Geodesy and Surveying in Python',
      url='https://github.com/GeoscienceAustralia/GeodePy',
      author='Geoscience Australia',
      author_email='geodesy@ga.gov.au',
      license='Apache License 2.0',
      packages=['geodepy'],
      install_requires=['numpy', 'scipy'])
