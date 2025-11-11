.. _features/gnss:

Sinex
=====

This module includes functions for dealing with sinex files.

Reading Sinex Files
----------------------

.. autofunction:: geodepy.gnss.list_sinex_blocks
.. autofunction:: geodepy.gnss.read_sinex_comments
.. autofunction:: geodepy.gnss.read_sinex_header_line
.. autofunction:: geodepy.gnss.read_sinex_custom
.. autofunction:: geodepy.gnss.read_sinex_estimate
.. autofunction:: geodepy.gnss.read_sinex_matrix
.. autofunction:: geodepy.gnss.read_sinex_sites
.. autofunction:: geodepy.gnss.read_disconts
.. autofunction:: geodepy.gnss.read_solution_epochs
.. autofunction:: geodepy.gnss.read_sinex_header_block
.. autofunction:: geodepy.gnss.read_sinex_file_reference_block
.. autofunction:: geodepy.gnss.read_sinex_input_acknowledgments_block
.. autofunction:: geodepy.gnss.read_sinex_solution_statistics_block
.. autofunction:: geodepy.gnss.read_sinex_site_receiver_block
.. autofunction:: geodepy.gnss.read_sinex_site_antenna_block
.. autofunction:: geodepy.gnss.read_sinex_site_gps_phase_center_block
.. autofunction:: geodepy.gnss.read_sinex_site_eccentricity_block
.. autofunction:: geodepy.gnss.read_sinex_site_id_block
.. autofunction:: geodepy.gnss.read_sinex_solution_epochs_block
.. autofunction:: geodepy.gnss.read_sinex_solution_estimate_block
.. autofunction:: geodepy.gnss.read_sinex_solution_apriori_block
.. autofunction:: geodepy.gnss.read_sinex_solution_matrix_estimate_block
.. autofunction:: geodepy.gnss.read_sinex_solution_matrix_apriori_block

Converting Sinex Data to DataFrames
------------------------------------

.. autofunction:: geodepy.gnss.matrix2dataframe_solution_matrix_estimate
.. autofunction:: geodepy.gnss.sinex2dataframe_solution_estimate
.. autofunction:: geodepy.gnss.sinex2dataframe_solution_apriori
.. autofunction:: geodepy.gnss.sinex2dataframe_solution_matrix_estimate

Writing Sinex Files
--------------------

.. autofunction:: geodepy.gnss.print_sinex_comments
.. autofunction:: geodepy.gnss.set_creation_time
.. autofunction:: geodepy.gnss.dataframe2sinex_solution_estimate
.. autofunction:: geodepy.gnss.dataframe2sinex_solution_apriori
.. autofunction:: geodepy.gnss.dataframe2matrix_snx_vcv
.. autofunction:: geodepy.gnss.dataframe2matrix_solution_matrix_estimate
.. autofunction:: geodepy.gnss.dataframe2sinex_solution_matrix_estimate
.. autofunction:: geodepy.gnss.writeSINEX

Specific Sinex Functions
--------------------------

.. autofunction:: geodepy.gnss.remove_stns_sinex
.. autofunction:: geodepy.gnss.remove_velocity_sinex
.. autofunction:: geodepy.gnss.remove_matrixzeros_sinex