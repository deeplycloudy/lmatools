from __future__ import absolute_import
from distutils.core import setup

setup(name='lmatools',
    version='0.1',
    description='Python tools for processing and visualizing Lightning Mapping Array data',
    author='Eric Bruning',
    author_email='eric.bruning@gmail.com',
    url='https://bitbucket.org/deeplycloudy/lmatools/',
    package_dir={'lmatools': ''}, # wouldn't be necessary if we reorganized to traditional package layout with lmatools at the same directory level as the setup.py script.
    packages=['lmatools', 'lmatools.flashsort', 'lmatools.flashsort.autosort'],
    )