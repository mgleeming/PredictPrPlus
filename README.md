# predictPrPlus

Introduction
------------------------------------------------------------------------------------

PredictPrPlus optimises the spatial charge distribution for multiply charged protein ions and estimates the apparent gas-phase basicity or acidity for protonated or deprotonated ions respectively. The optimal charge configuration depends upon the 'intrinsic' basicity/acidity of the residues in a protein as well as the distance between charges. Optimisation of the configuration is performed using the 'pseudo-random walk' algorithm described by Schnier *et al* (*J. Am. Soc. Mass Spectrom.* 1995, 6, 1086-1097).

Arbitrary peptides can be modelled as linear strings of residues which provide reasonable approximations of protein ions at high charge states. Alternatively, gas-phase basicity/acidity values can be determined for more complex structures by providing a PDB structure file.


Installation
------------------------------------------------------------------------------------
	sudo python setup.py build_ext --inplace install

predictPrPlus relies on the following libraries:
- numpy
- cython
- qt4
- PyQt4
- pyqtgraph
- prody

Parameters and usage
------------------------------------------------------------------------------------

#### Input structures
PredictPrPlus can optimise charge distributions for proteins as either linear strings or three-dimensional structures. If optimisation of a linear string is required, enter the protein primary sequence (one letter amino acid code) into the 'Protein Sequence' text box. Alternatively, for 3D structures, specify the desired protein structure by loading a PDB file. In this case, it is important that atom and residue labels in the PDB file conform to PDB naming conventions.

In some cases, it may be desirable to include a fixed charge site in the simulation (i.e. a site that remains charged in all cases). For example, this can may represent the effects of a coordinating metal ion. When modelling a protein as a linear string, simply add an 'X' at the desired fixed charge site. If a PDB file is specified, specify the desired cartesian coordinates of the fixed charge site using the 'Add Fixed Charge' tool.

#### Optimisation parameters:

| Parameter 								| Description																											     |
| -------------------------	| -------------------------------------------------------------------- |
| Dielectric Constant       | Dielectric shielding of electrostatic effects between charges   		 |
| Residue Spacing 					| Length of one amino acid residue in Angstroms (for linear only) 		 |
| Energy Cutoff							| Store charge configurations within this energy of the global minimum |
| Optimisation steps				| Number of charge movements performed per iteration 									 |
| Iterations								| Number of iterations to perform 																		 |
| Ion polarity 							| Simulate cations (gas-phase basicity) or anions (gas-phase acidity)	 |
| Charge state range				| Range of charge states to optimise 																	 |

The intrinsic gas-phase basicity or gas-phase acidity values used in the optimisation process can be defined by clicking 'Edit GB/GA intrinsic values'. For amino acid residues, values indicate the GB/GA of the residue side chain. If the field is blank, no ionisable side chain is defined for that residue. The GB/GA values for backbone amide bonds are the same for all residues and can be altered by changing the values in the 'amide' field. The GB/GA values for protein N-terminal amine and C-terminal carboxylate groups can be specified using the 'N-term' and 'C-term' fields respectively.


License
------------------------------------------------------------------------------------

    Copyright (c) 2016
    Michael G. Leeming. William A. Donald, Muhammad Bin Zenaidee.

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the author may not be used to endorse or promote products
       derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
    IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
    NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
    THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
