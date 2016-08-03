from __future__ import absolute_import
<<<<<<< HEAD
from LMAtools.multiples_nc import make_plot
=======
from LMAtools.vis.multiples_nc import make_plot
>>>>>>> 877fb32c26c34947a91aa3e27f815967bc9a17f2

n_cols=2

make_plot('LMA_20090724_010000_240_flash_extent.nc', 'flash_extent', n_cols=n_cols)
make_plot('LMA_20090724_010000_240_flash_init.nc', 'flash_initiation', n_cols=n_cols)
make_plot('LMA_20090724_010000_240_source.nc', 'lma_source', n_cols=n_cols)
make_plot('LMA_20090724_010000_240_footprint.nc', 'flash_footprint', n_cols=n_cols)