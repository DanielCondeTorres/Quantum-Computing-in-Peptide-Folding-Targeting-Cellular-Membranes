# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""A class defining a Miyazawa-Jernigan interaction between beads of a peptide."""

import numpy as np

from .interaction import Interaction
from ..residue_validator import _validate_residue_sequence
from ..data_loaders.energy_matrix_loader import (
    _load_energy_matrix_file,
)

# pylint: disable=too-few-public-methods


class MiyazawaJerniganInteraction(Interaction):
    """A class defining a Miyazawa-Jernigan interaction between beads of a peptide.
    Details of this model can be found in Miyazawa, S. and Jernigan, R. L. J. Mol. Biol.256,
    623–644 (1996), Table 3."""

    def calculate_energy_matrix(self, residue_sequence: str) -> np.ndarray: # ,medio_w_1: float, medio_cs_1: float, exchange_solvent_solvent: float) -> np.ndarray:
        """
        Calculates an energy matrix for a Miyazawa-Jernigan interaction based on the
        Miyazawa-Jernigan potential file.

        Args:
            residue_sequence: A string that contains characters defining residues for
                            a chain of proteins.

        Returns:
            Numpy array of pair energies for amino acids.
        """
        chain_len = len(residue_sequence)
        _validate_residue_sequence(residue_sequence)
        mj_interaction, list_aa = _load_energy_matrix_file()
        pair_energies = np.zeros((chain_len + 1, 2, chain_len + 1, 2))
        for i in range(1, chain_len + 1):
            for j in range(i + 1, chain_len + 1):
                aa_i = list_aa.index(residue_sequence[i - 1])
                aa_j = list_aa.index(residue_sequence[j - 1])
                pair_energies[i, 0, j, 0] = mj_interaction[(min(aa_i, aa_j)), max(aa_i, aa_j)]
        return pair_energies



    def calculate_media_energy_matrix(self, residue_sequence: str,phase_w_1: float, phase_cs_1: float, phase_w_2:float, phase_cs_2:float, exchange_solvent_solvent: float, correction_mj: bool):

        
        '''       
                This function modifies the Miyazawa-Jernigan interaction to take into account in which medium each bead is located.
                Details of the model can be found in Leonhard, K., Prausnitz, J.M. and Radke, C.J. (2004), 
                Solvent–amino acid interaction energies in three-dimensional-lattice Monte Carlo simulations of a model 27-mer protein: Folding thermodynamics and kinetics. 
                Protein Science, 13: 358-369. https://doi.org/10.1110/ps.03198204

                Inputs: 

                * residue_sequence: amino acid sequence
                * ws_phase_1: average interaction of an amino acid with solvent in phase 1
                * ws_phase_2: average interaction of an amino acid with solvent in phase 2
                * cs_phase_1: model solvent in phase 1
                * cs_phase_2: model solvent in phase 2
                * exchange_solvent: energy exchange between interfaces
                * correction_mj: indicates whether we decide to use the Miyazawa-Jerningan potential upgrade or not.

                Output:

                * solvent_phase_1: Numpy array of pair energies for amino acids if is in phase 1
                * solvent_phase_2: Numpy array of pair energies for amino acids if is in phase 2
                * exchange_solvent_solvent: Energy exchange between interfaces
        '''
        
        # It checks if the correction_mj flag is set to True. If so, it will proceed with the calculations: corrects Miyazawa-Jerningan potential.
        if correction_mj == True:
            
            # nitialize Sequence Length & Validate
            chain_len = len(residue_sequence)
            _validate_residue_sequence(residue_sequence)

            # Load Energy Matrix:
            mj_interaction, list_aa = _load_energy_matrix_file()
            
            # Initialises the two possible solvent arrays to be selected according to the phase we are in.
            solvent_phase_1 = np.zeros((chain_len + 1, 2, chain_len + 1, 2))
            solvent_phase_2 = np.zeros((chain_len + 1, 2, chain_len + 1, 2))
            # Compute Interactions:
            for i in range(1, chain_len + 1):
                aa_i = list_aa.index(residue_sequence[i - 1])
                solvent_phase_1[i, 0, i, 0]=1/2*(1-phase_cs_1)*mj_interaction[min(aa_i, aa_i),max(aa_i, aa_i)]+phase_w_1
                solvent_phase_2[i, 0, i, 0]=1/2*(1-phase_cs_2)*mj_interaction[min(aa_i, aa_i),max(aa_i, aa_i)]+phase_w_2
                
                # A nested loop then computes the sum of energies (ejj) for the amino acid with all other amino acids in the chain: sum_{j}^{N}=ejj
                ejj = 0
                for j in range(1, chain_len + 1):#maybe the 20 aa
                    aa_j = list_aa.index(residue_sequence[j - 1])
                    ejj += ejj + mj_interaction[min(aa_j, aa_j),max(aa_j, aa_j)]
                
                # Result for each phase: Updated 
                #change chain_len with 20 aa
                solvent_phase_1[i, 0, i, 0]= solvent_phase_1[i, 0, i, 0] + phase_cs_1/(2*chain_len)*ejj
                solvent_phase_2[i, 0, i, 0]= solvent_phase_2[i, 0, i, 0] + phase_cs_2/(2*chain_len)*ejj
                print('SELF 1',solvent_phase_1, solvent_phase_2)
            return solvent_phase_1,solvent_phase_2,exchange_solvent_solvent
            
        # If not, it will return None values: This implies working with the original Miyazawa-Jerningan potential.
        else:
            return None, None, None


