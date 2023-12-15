import numpy as np 
from qiskit.opflow import OperatorBase, PauliOp, PauliSumOp

from typing import Union
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from ..peptide.pauli_ops_builder import _build_full_identity

class Auxiliar_mj_modification:
    
    def calculate_pairs_mj_modification(i: int, j: int , pair_energies_phase_1: np.ndarray , pair_energies_phase_2: np.ndarray, sign_function, weight_interface: float, energy_exchange: float, displacement_of_the_interface_plane: float):
        print('PARES: ',i,j)
        energy_phase_1_bead_1 = pair_energies_phase_1[i][0][i][0]
        energy_phase_1_bead_2 = pair_energies_phase_1[j][0][j][0]
        energy_phase_2_bead_1 = pair_energies_phase_2[i][0][i][0]
        energy_phase_2_bead_2 = pair_energies_phase_2[j][0][j][0]
  
        # Adjust the energies based on the positions of beads relative to the interface
        if i > 0:
            lower_bead = sign_function[i-2]
            upper_bead = sign_function[j-2]
  
            # Calculate phase energy contributions for beads near the interfacei
  
            '''
            Here we are creating terms such that: if the sign function is positive (negative), and the value of the interface is positive (negative) we get a 0 to say that the bead is in phase 2,
            however, we will have a 1 for phase 1. As we have the code the positive (negative) interface corresponds to a hydrophobic (hydrophilic) phase 1 and a hydrophilic (hydrophobic) phase 2,
            thus leaving phase 1 up-hydrophobic (down-hydrophilic):
            '''
  
            energy_phase_bead_lower = -1*(lower_bead - weight_interface/abs(weight_interface) *_build_full_identity(lower_bead.num_qubits))*energy_phase_2_bead_1 + (lower_bead + weight_interface/abs(weight_interface)*_build_full_identity(lower_bead.num_qubits))*energy_phase_1_bead_1
  
            energy_phase_bead_upper = -1*(upper_bead - weight_interface/abs(weight_interface)*_build_full_identity(upper_bead.num_qubits))*energy_phase_2_bead_2 + (upper_bead + weight_interface/abs(weight_interface)*_build_full_identity(upper_bead.num_qubits))*energy_phase_1_bead_2
  
            # Combine energy contributions from different beads
            energy_phase_term = 0.5*(energy_phase_bead_lower+energy_phase_bead_lower) - energy_exchange*_build_full_identity(lower_bead.num_qubits)
  
        elif i == 0:
            upper_bead = sign_function[j-2]
            Identity = _build_full_identity(upper_bead.num_qubits)
            if displacement_of_the_interface_plane != 0:
                IP = displacement_of_the_interface_plane/abs(displacement_of_the_interface_plane)
            else:
                IP = 0
            # Calculate phase energy contribution for the first bead
            '''
            Here, if IP = -1 (1), our first bead is below (on top of) the plane if the value of weight interface is negative (positive) we get a 0 in phase 2, so our bead have to be in phase 1
  
            '''
            energy_phase_first_bead = -1*(IP- weight_interface/abs(weight_interface))*energy_phase_2_bead_1 + (IP + weight_interface/abs(weight_interface))*energy_phase_1_bead_1
  
            energy_phase_bead_upper = -1*(upper_bead - weight_interface/abs(weight_interface)*_build_full_identity(upper_bead.num_qubits))*energy_phase_2_bead_2 + (upper_bead + weight_interface/abs(weight_interface)*_build_full_identity(upper_bead.num_qubits))*energy_phase_1_bead_2
  
            # Combine energy contributions
            print('energy first bead',energy_phase_first_bead)
  
            energy_phase_term = 0.5*(energy_phase_first_bead*Identity+energy_phase_bead_upper) - energy_exchange*Identity
        return energy_phase_term
