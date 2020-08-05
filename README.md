# opiateFOSconnectivity

Citation:
Brynildsen JK, Mace KD, Cornblath EJ, Weidler C, Pasqualetti F, Bassett DS, Blendy JA (2020) Gene coexpression patterns predict opiate-induced brain-state transitions. Proceedings of the National Academy of Sciences (https://www.pnas.org/content/early/2020/07/20/2003601117)

Data files:
.csv files contain the fold change in FOS expression data for each condition described in the manuscript (Fig. 2).

To obtain gene expression data from the Allen Institute, please see https://github.com/ejcorn/mouse_abi_tool

Code:
FOS_corr_network.R - used to generate FOS correlation networks and compute graph theory metrics (Fig. 3)
Gene_Coexpression.R - used to analyze normalized gene expression data from the Allen Institute (Fig. 4)
Emin_diff_node.m, min_eng_cont.m - used to compute change in minimum control energy (Fig. 5)
