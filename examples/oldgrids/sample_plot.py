from __future__ import absolute_import
from LMAtools.vis.multiples_nc import make_plot

n_cols=2

make_plot('LMA_20090724_010000_240_flash_extent.nc', 'flash_extent', n_cols=n_cols)
make_plot('LMA_20090724_010000_240_flash_init.nc', 'flash_initiation', n_cols=n_cols)
make_plot('LMA_20090724_010000_240_source.nc', 'lma_source', n_cols=n_cols)
make_plot('LMA_20090724_010000_240_footprint.nc', 'flash_footprint', n_cols=n_cols)
