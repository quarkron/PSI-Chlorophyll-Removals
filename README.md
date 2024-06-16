This requires an installation of UCSF Chimera 2 (not Chimera X!)

Steps:
1. Run chl_extractor_v2.py in the Python IDLE interface of UCSF Chimera. This should write a 'chls_1B0.csv' file in the directory containing chlorophyll positions.
2. Run rate network_arrows.py in Chimera. This creates a visualization of the transfer rate network.
3. Open node_metrics_calculator.py using your Python installation. Set create_data to True. This creates 'PSI_chls.csv' file containing the node metrics per chlorophyll. Alternatively, run data_gen.ipynb.
4. Run cla_characteristics_v2.py in Chimera to visualize the nodes with their metric values.
5. Run qeffBC_datagen.ipynb, neffBC_datagen.ipynb, qeffRandom_datagen.ipynb, neffRandom_datagen.ipynb. The neff files will take a while (may take days depending on the number of runs). Modify the files accordingly to parallelize the computations.
6. Run fig_gen2.ipynb to generate the qeff, neff, robustness plots.

Notes:
The graph_constructor.py and node_metrics_calculator.py contain the bulk of the source code. Algorithms and graph-theoretic tools can be found in graph_constructor. 

If you have any questions, please don't hesitate to email me at the following email addresses:
rmacda@nip.edu.upd.ph
racda@illinois.edu 
acdaronmichael@gmail.com
