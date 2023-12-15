# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Builds qubit operators for all Hamiltonian terms in the protein folding problem."""
from typing import Union
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np
from qiskit.opflow import OperatorBase, PauliOp, PauliSumOp
from .peptide.chains.main_chain import BaseChain
from .bead_contacts.contact_map import ContactMap
from .bead_distances.distance_map import DistanceMap
from .exceptions.invalid_side_chain_exception import (InvalidSideChainException,)
from .exceptions.invalid_size_exception import InvalidSizeException
from .penalty_parameters import PenaltyParameters
from .peptide.pauli_ops_builder import _build_full_identity
from .qubit_utils.qubit_fixing import _fix_qubits
from .peptide.beads.base_bead import BaseBead
from .peptide.peptide import Peptide

from .interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction
from .interface_bead_interaction.interface_parameters import  interface_parameters
from .interface_bead_interaction.auxiliar_mj_modification_hamiltonian import Auxiliar_mj_modification
# pylint: disable=too-few-public-methods

class QubitOpBuilder:
    """Builds qubit operators for all Hamiltonian terms in the protein folding problem."""
    def __init__(self, peptide: Peptide, pair_energies: np.ndarray, penalty_parameters: PenaltyParameters, pair_energies_modified_mj: np.ndarray , cs_phase_1: float, cs_phase_2: float, w_phase_1:float, w_phase_2: float, one_solvent_parameter: bool):
        """Builds qubit operators for all Hamiltonian terms in the protein folding problem.

        Args:
            peptide: A Peptide object that includes all information about a protein.
            pair_energies: Numpy array of pair energies for amino acids.
            penalty_parameters: A PenaltyParameters object storing the values of all penalty parameters.
        """
        
        # Energy modification MJ potential parameters
        self._pair_energies_phase_1, self._pair_energies_phase_2, self._exchange = pair_energies_modified_mj 
        self._cs_phase_1 = cs_phase_1; self._cs_phase_2 = cs_phase_2
        self._w_phase_1 = w_phase_1; self._w_phase_2 = w_phase_2
        self._one_solvent_parameter = one_solvent_parameter
        self._interface_value = interface_parameters().calculate_weight_interface (self._cs_phase_1, self._cs_phase_2, self._w_phase_1 ,self._w_phase_2)
        
        self._peptide = peptide
        self._pair_energies = pair_energies
        self._penalty_parameters = penalty_parameters
        self._contact_map = ContactMap(peptide)
        self._distance_map = DistanceMap(peptide)
        _side_chain_hot_vector = self._peptide.get_side_chain_hot_vector()
        self._has_side_chain_second_bead = _side_chain_hot_vector[1] if len(_side_chain_hot_vector) > 1 else False
    
    def build_qubit_op(self) -> Union[PauliSumOp, PauliOp]:
        """
        Builds a qubit operator for a total Hamiltonian for a protein folding problem. It includes
        8 terms responsible for chirality, geometry and nearest neighbors interactions.

        Returns:
            A total Hamiltonian for the protein folding problem.

        Raises:
            InvalidSizeException: if chains of invalid/incompatible sizes provided.
            InvalidSideChainException: if side chains on forbidden indices provided.
        """
        side_chain = self._peptide.get_side_chain_hot_vector()
        main_chain_len = len(self._peptide.get_main_chain)

        if len(side_chain) != main_chain_len: raise InvalidSizeException("side_chain_lens size not equal main_chain_len")
        if side_chain[0] == 1 or side_chain[-1] == 1:
            raise InvalidSideChainException(
                "First and last main beads are not allowed to have a side chain. Nonempty "
                "residue provided for an invalid side chain."
            )

        num_qubits = 4 * pow(main_chain_len - 1, 2)
        full_id = _build_full_identity(num_qubits) # build a full identity operator of dimension num_qubits I**num_qubits

        h_chiral = self._create_h_chiral()
        if h_chiral != 0: h_chiral = full_id ^ h_chiral
        
        h_back = self._create_h_back()
        if h_back != 0: h_back = full_id ^ h_back

        h_scsc = self._create_h_scsc() if self._penalty_parameters.penalty_1 else 0
        h_bbbb = self._create_h_bbbb() if self._penalty_parameters.penalty_1 else 0

        h_short = self._create_h_short()
        if h_short != 0: h_short = full_id ^ h_short

        h_bbsc, h_scbb = self._create_h_bbsc_and_h_scbb() if self._penalty_parameters.penalty_1 else (0, 0)

        h_interface, sign_function= self._create_h_interface()
        if h_interface != 0: h_interface = full_id ^ h_interface
        
        
        h_bbbb_media = self._create_h_bbbb_media_dependent() if self._penalty_parameters.penalty_1 else 0
        

        h_total = h_chiral + h_back + h_short + h_bbbb + h_bbsc + h_scbb + h_scsc + h_interface + h_bbbb_media
        
        return h_total.reduce()

    def _create_turn_operators(self, lower_bead: BaseBead, upper_bead: BaseBead) -> OperatorBase:
        """
        Creates a qubit operator for consecutive turns.

        Args:
            lower_bead: A bead with a smaller index in the chain.
            upper_bead: A bead with a bigger index in the chain.

        Returns:
            A qubit operator for consecutive turns.
        """
        lower_bead_indic_0, lower_bead_indic_1, lower_bead_indic_2, lower_bead_indic_3 = lower_bead.indicator_functions
        upper_bead_indic_0, upper_bead_indic_1, upper_bead_indic_2, upper_bead_indic_3 = upper_bead.indicator_functions

        turns_operator = _fix_qubits(
            lower_bead_indic_0 @ upper_bead_indic_0
            + lower_bead_indic_1 @ upper_bead_indic_1
            + lower_bead_indic_2 @ upper_bead_indic_2
            + lower_bead_indic_3 @ upper_bead_indic_3,
            self._has_side_chain_second_bead,
        )
        return turns_operator

    def _create_turn_operators_axis(self, lower_bead: BaseBead, upper_bead: BaseBead, axis:int )-> OperatorBase:
        """
        Crea un operador de qubits sobre un unico eje.

        Args:
            lower_bead: Bead inferior, empezaremos usando la 0
            upper_bead: Bead superiores, se calculara para cada bead
            axis: eje sobre el que esstamos trabajando

        Return:
            Operador de qubits para un eje especifico
        """
        lower_bead_indic_0, lower_bead_indic_1, lower_bead_indic_2, lower_bead_indic_3 = lower_bead.indicator_functions
        upper_bead_indic_0, upper_bead_indic_1, upper_bead_indic_2, upper_bead_indic_3 = upper_bead.indicator_functions
        diccionario_eje_seleccionado = {
            0: (lower_bead_indic_0, upper_bead_indic_0), #el primer valor es el eje que seleccionamos
            1: (lower_bead_indic_1, upper_bead_indic_1),
            2: (lower_bead_indic_2, upper_bead_indic_2),
            3: (lower_bead_indic_3, upper_bead_indic_3)
            }
        par_seleccionado = diccionario_eje_seleccionado.get(axis)
        turn_operator_axis=_fix_qubits(par_seleccionado[0] @ par_seleccionado[1])
        return turn_operator_axis

    def _create_h_back(self) -> Union[PauliSumOp, PauliOp]:
        """
        Creates Hamiltonian that imposes the geometrical constraint wherein consecutive turns along
        the same axis are penalized by a factor, penalty_back. Note, that the first two turns are
        omitted (fixed in optimization) due to symmetry degeneracy.

        Returns:
            Contribution to Hamiltonian in symbolic notation that penalizes consecutive turns
            along the same axis.
        """

        main_chain = self._peptide.get_main_chain
        penalty_back = self._penalty_parameters.penalty_back
        h_back = 0
        for i in range(len(main_chain) - 2):
            h_back += penalty_back * self._create_turn_operators(main_chain[i], main_chain[i + 1])
        h_back = _fix_qubits(h_back, self._has_side_chain_second_bead)
        return h_back

    def _create_h_chiral(self) -> Union[PauliSumOp, PauliOp]:
        """
        Creates a penalty/constrain term to the total Hamiltonian that imposes that all the position
        of all side chain beads impose the right chirality. Note that the position of the side chain
        bead at a location (i) is determined by the turn indicators at i - 1 and i. In the absence
        of side chains, this function returns a value of 0.

        Returns:
            Hamiltonian term that imposes the right chirality.
        """

        main_chain = self._peptide.get_main_chain
        main_chain_len = len(main_chain)
        h_chiral = 0
        # 2 stands for 2 qubits per turn, another 2 stands for main and side qubit register
        full_id = _build_full_identity(2 * 2 * (main_chain_len - 1))
        for i in range(1, len(main_chain) + 1):
            upper_main_bead = main_chain[i - 1]
            if upper_main_bead.side_chain is None: continue
            upper_side_bead = upper_main_bead.side_chain[0]
            lower_main_bead = main_chain[i - 2]

            lower_main_bead_indic_0, lower_main_bead_indic_1, lower_main_bead_indic_2, lower_main_bead_indic_3 = lower_main_bead.indicator_functions
            upper_main_bead_indic_0, upper_main_bead_indic_1, upper_main_bead_indic_2, upper_main_bead_indic_3 = upper_main_bead.indicator_functions
            upper_side_bead_indic_0, upper_side_bead_indic_1, upper_side_bead_indic_2, upper_side_bead_indic_3 = upper_side_bead.indicator_functions
            turn_coeff = int((1 - (-1) ** i) / 2)
            h_chiral += self._build_chiral_term(full_id, lower_main_bead_indic_1, lower_main_bead_indic_2, lower_main_bead_indic_3, turn_coeff,
                upper_main_bead_indic_1, upper_main_bead_indic_2, upper_main_bead_indic_3, upper_side_bead_indic_0,)
            h_chiral += self._build_chiral_term(full_id, lower_main_bead_indic_0, lower_main_bead_indic_3, lower_main_bead_indic_2, turn_coeff,
                upper_main_bead_indic_0, upper_main_bead_indic_3, upper_main_bead_indic_2, upper_side_bead_indic_1,)
            h_chiral += self._build_chiral_term(full_id, lower_main_bead_indic_0, lower_main_bead_indic_1, lower_main_bead_indic_3, turn_coeff,
                upper_main_bead_indic_0, upper_main_bead_indic_1, upper_main_bead_indic_3, upper_side_bead_indic_2,)
            h_chiral += self._build_chiral_term(full_id, lower_main_bead_indic_0, lower_main_bead_indic_2, lower_main_bead_indic_1, turn_coeff,
                upper_main_bead_indic_0, upper_main_bead_indic_2, upper_main_bead_indic_1, upper_side_bead_indic_3,)
            h_chiral = _fix_qubits(h_chiral, self._has_side_chain_second_bead)
        return h_chiral

    def _build_chiral_term(self, full_id, lower_main_bead_indic_b, lower_main_bead_indic_c, lower_main_bead_indic_d,turn_coeff,
        upper_main_bead_indic_b, upper_main_bead_indic_c, upper_main_bead_indic_d, upper_side_bead_indic_a,):
        
        lmbb=lower_main_bead_indic_b; lmbc=lower_main_bead_indic_c; lmbd=lower_main_bead_indic_d
        umbb=upper_main_bead_indic_b; umbc=upper_main_bead_indic_c; umbd=upper_main_bead_indic_d; usba=upper_side_bead_indic_a
        
        return (self._penalty_parameters.penalty_chiral * (full_id - upper_side_bead_indic_a)
            @ ((1 - turn_coeff) * (lmbb @ umbc + lmbc @ umbd + lmbd @ umbb)
                  + turn_coeff  * (lmbc @ umbb + lmbd @ umbc + lmbb @ umbd)))
    
    def _create_h_bbbb(self) -> Union[PauliSumOp, PauliOp]:
        """
        Creates Hamiltonian term corresponding to a 1st neighbor interaction between main/backbone (BB) beads.

        Returns:
            Hamiltonian term corresponding to a 1st neighbor interaction between main/backbone (BB) beads.
        """
        penalty_1 = self._penalty_parameters.penalty_1
        h_bbbb = 0
        main_chain_len = len(self._peptide.get_main_chain)
        
        # Use this Hamiltoninan if we want the original MJ potential in the interface, if only one phase is selected depends on the correction term
        if self._pair_energies_phase_1 is None or self._interface_value == 0:# the second option corresponds to self.pair_energies_phase_1 =self.pair_energies_phase_2
            if self._pair_energies_phase_1 is None:
                _pair_energies_phase_1 = 0*self._pair_energies

            else:
                _pair_energies_phase_1 = self._pair_energies_phase_1

            for i in range(1, main_chain_len-3):
                for j in range(i + 5, main_chain_len + 1):
                    if (j - i) % 2 == 0: continue
                    energy_phase_1_bead_1 = _pair_energies_phase_1[i][0][i][0]
                    energy_phase_1_bead_2 = _pair_energies_phase_1[j][0][j][0]
                    modification_mj =  energy_phase_1_bead_1 + energy_phase_1_bead_2
                    h_bbbb += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.first_neighbor(self._peptide, i, 0, j, 0, penalty_1, self._pair_energies, modification_mj))
                    

                    try:
                        energy_phase_1_bead_1 = _pair_energies_phase_1[i-1][0][i-1][0]
                        energy_phase_1_bead_2 = _pair_energies_phase_1[j][0][j][0]
                        modification_mj =  energy_phase_1_bead_1 + energy_phase_1_bead_2
                        h_bbbb += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide,i - 1,0,j,0,penalty_1,self._pair_energies, modification_mj))
                    except (IndexError, KeyError): pass
                    
                    try:
                        energy_phase_1_bead_1 = _pair_energies_phase_1[i+1][0][i+1][0]
                        energy_phase_1_bead_2 = _pair_energies_phase_1[j][0][j][0]
                        modification_mj =  energy_phase_1_bead_1 + energy_phase_1_bead_2
                        h_bbbb += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide,i + 1,0,j,0,penalty_1,self._pair_energies, modification_mj))
                    except (IndexError, KeyError): pass
                    
                    try:
                        energy_phase_1_bead_1 = _pair_energies_phase_1[i][0][i][0]
                        energy_phase_1_bead_2 = _pair_energies_phase_1[j-1][0][j-1][0]
                        modification_mj =  energy_phase_1_bead_1 + energy_phase_1_bead_2
                        h_bbbb += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide,i,0,j - 1,0,penalty_1,self._pair_energies, modification_mj))
                    except (IndexError, KeyError): pass
                    
                    try:
                        energy_phase_1_bead_1 = _pair_energies_phase_1[i][0][i][0]
                        energy_phase_1_bead_2 = _pair_energies_phase_1[j+1][0][j+1][0]
                        modification_mj =  energy_phase_1_bead_1 + energy_phase_1_bead_2
                        h_bbbb += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide,i,0,j + 1,0,penalty_1,self._pair_energies,modification_mj))
                    except (IndexError, KeyError): pass
                    h_bbbb = _fix_qubits(h_bbbb, self._has_side_chain_second_bead)
        else:
            h_bbbb=0
        return h_bbbb
    
    def _create_h_interface(self) -> Union[PauliSumOp, PauliOp]: 

        """
            Hamiltonian term corresponding to the produced by the interface taking into account the hydrophobicity of the different amino acids for the different media.    
        """
        # Access the dictionary containing information on the hydrophobicity of the different amino acids.
        amino_info =  interface_parameters.AMINO

        # We call the class containing the amino acid sequence.
        main_chain_residue_sequence = self._peptide.get_main_chain.main_chain_residue_sequence # Peptide sequence.
        main_chain_len = len(self._peptide.get_main_chain) # Peptide sequence's length. 
        
        # Empty list to store the hydrophobicity values of our aa sequence.
        hydrophobicity_value=[] 
        
        
        for target_letter in main_chain_residue_sequence: # Walk through the amino acid sequence.
            # Look up the letter in the dictionary and get hydrophobicity if there is a match.
            for key,value in amino_info.items(): 
                if value['label'] == target_letter: # 'label' stores the letters in the amino_info dictionary, we look up the entry for the desired letter
                    hydrophobicity_value.append(value['hydrophobicity'])  #We obtain the value of hydrophobicity for that letter.

        full_id = _build_full_identity(2 * 2 * (main_chain_len - 1))       
        penalty_1 = self._penalty_parameters.penalty_1
        h_interface = 0
        displacementMap = self._distance_map.displacement_map(interface_parameters.axis)
        

        
        bead_initial_distance = self._peptide.get_main_chain[0]

        distanceMap = self._distance_map.distance_map
    
        beads = list(displacementMap.keys()) # We call the information contained in all the beads of our peptide.
        first_bead=beads[0] # As the first bead is fixed, we are interested in it to calculate the distances with respect to it.
        
        #Cambiar por la identidad, mejor que hacer el bucle emplear el full_id
        i = -1
        for bead, value in displacementMap[first_bead].items(): 
            if bead: 
                i += 1
                if (i==1):
                    one = value.reduce()
                    break
        print("one = ", one)
        print("Ide = ",-1*full_id ) #Creo que esto es mejor, no tenemos que hacer bucle ni depender de obligar a que sea en el axis 1.
        one = -1*full_id # Para tenerlo bien ya
        print("\n\n Now we are going to print the displacementMap that will be used in the interface hamiltonian: \n\n",displacementMap[first_bead].items())
        TotalDisplacements = [] # Empty list to store the displacements of each bead with respect to the initial one.
        Sign_function_value = [] # Empty list to store the values (aproximated) of the Sing function polynomial aproximation. 
        
        i = -1 # Start counter
        w = 0
        for bead, value in displacementMap[first_bead].items(): # We get the distance to the initial bead
            if bead: # If bead exists  
                i += 1
                print('Same: ',(value.reduce(),'that',interface_parameters.displacement_of_the_interface_plane*one))
                TotalDisplacements.append((value.reduce()-interface_parameters.displacement_of_the_interface_plane*one).reduce()) # We obtain the displacement of each bead with respect to the plane, not with respect to the initial bead, thanks to this displacement of the plane.  
                
                # Approximation of the sign function to a polynomial expression.
                TD = TotalDisplacements[i]
                TD2 = (TD@TD).reduce()
                TD3 = (TD@TD2).reduce()
                TD4 = (TD2@TD2).reduce()
                TD5 = (TD2@TD3).reduce()
                TD6 = (TD3@TD3).reduce()
                TD7 = (TD3@TD4).reduce()
                # Coefficents
                coef_7 = interface_parameters.coefficients_of_chebysev['coef_7']; coef_6 = interface_parameters.coefficients_of_chebysev['coef_6']   
                coef_5 = interface_parameters.coefficients_of_chebysev['coef_5']; coef_4 = interface_parameters.coefficients_of_chebysev['coef_4']
                coef_3 = interface_parameters.coefficients_of_chebysev['coef_3']; coef_2 = interface_parameters.coefficients_of_chebysev['coef_2']
                coef_1 = interface_parameters.coefficients_of_chebysev['coef_1']; coef_0 = interface_parameters.coefficients_of_chebysev['coef_0']
                #Result
                sign_function_approximation = coef_7*TD7+coef_6*TD6+coef_5*TD5+coef_4*TD4+coef_3*TD3+coef_2*TD2+coef_1*TD+coef_0*full_id
                Sign_function_value.append((sign_function_approximation.reduce())) # Safe sing function approximation for each aminoacid 
                
            
                # Obtain the h_interface contribution for each aminoacid. This function will be minimize
            
                ''' CHOOSE OPTION 1 OR 2 FOR ONE SOLVENT: These options do not affect the calculation of two phases '''
                
               
               #Option 1: We do not take into account the hydrophobicity parameter of the solvent in a phase. Will give a 0, so original MJ is recover
                if  self._one_solvent_parameter == False:
                    h_interface +=  self._interface_value * sign_function_approximation * hydrophobicity_value[i]
                    #Finish option 1
                #Option 2 We take into account the hydrophobicity parameter of the solvent in a phase:
                elif  self._one_solvent_parameter == True:
                    if self._interface_value != 0:
                        h_interface +=  self._interface_value * sign_function_approximation * hydrophobicity_value[i]
                    else: #In this option csi =csj
                        # Negative is more stable, so if aa is hydrophobic ( > 0) will interact with the interphace, remember * cs < 0: model hydrophobic solvents (oil), so that should be negative
                        # and also  w < 0:  aminoacids interact on average attractively with the solvent
                        w += 1
                        upper_bead =  self._peptide.get_main_chain[w]
                        print('distancia: ',distanceMap[bead_initial_distance][upper_bead], 'bead: ',upper_bead, bead_initial_distance)
                        # no tocar este 0.1 para que sea siempre igual
                        h_interface +=  hydrophobicity_value[i]* self._cs_phase_1*(distanceMap[bead_initial_distance][upper_bead]**2)#Hidrofobicidad marcada por self.cs_phase_1 y tambien si es alta o baja
                #Finish option 2
        return _fix_qubits(h_interface, self._has_side_chain_second_bead), Sign_function_value
    
    

    def _create_h_bbbb_media_dependent(self):  # -> Union[PauliSumOp, PauliOp]:
        """
        This function creates a media-dependent Hamiltonian term for a specific peptide system.
        We need an interface, to see when we need one matrix energy solvent (phase_1) or (phase_2)
        """

        main_chain_len = len(self._peptide.get_main_chain)  # Get the length of the main chain of the peptide
        h_interface, sign_function = self._create_h_interface()  # Get interface Hamiltonian and sign function
        penalty_1 = self._penalty_parameters.penalty_1  # Extract penalty parameter
        h_bbbb_media = 0  # Initialize Hamiltonian term for media
        main_chain_len = len(self._peptide.get_main_chain)  # (Note: Redundant line. The length was calculated above)
        energy_exchange = self._exchange  # Extract energy exchange parameter
        weight_interface = self._interface_value  # Get the weight of the interface
        plane_displacement = interface_parameters.displacement_of_the_interface_plane
        # Check if there are phase 1 energies and the interface weight is not zero
        if self._pair_energies_phase_1 is not None and weight_interface != 0:
            print('SELF :', self._pair_energies_phase_1, self._pair_energies_phase_2)
            for i in range(1, main_chain_len-3):
                for j in range(i + 5, main_chain_len + 1):
                    # Consider only pairs (i, j) with odd j-i difference
                    if (j-i) % 2 == 0: continue
                    energy_phase_term = Auxiliar_mj_modification.calculate_pairs_mj_modification(i, j, self._pair_energies_phase_1, self._pair_energies_phase_2, sign_function, weight_interface, energy_exchange, plane_displacement)
                    

                    # Add contributions to the main Hamiltonian term
                    h_bbbb_media += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.first_neighbor(self._peptide, i, 0, j, 0, penalty_1, self._pair_energies, energy_phase_term))
                
                    # Add contributions from neighboring beads if they exist
                    try:
                        energy_phase_term = Auxiliar_mj_modification.calculate_pairs_mj_modification(i-1, j, self._pair_energies_phase_1, self._pair_energies_phase_2, sign_function, weight_interface, energy_exchange, plane_displacement)
                        h_bbbb_media += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide, i - 1, 0, j, 0, penalty_1, self._pair_energies, energy_phase_term))
                    except (IndexError, KeyError): pass
                
                    try:
                        energy_phase_term = Auxiliar_mj_modification.calculate_pairs_mj_modification(i+1, j, self._pair_energies_phase_1, self._pair_energies_phase_2, sign_function, weight_interface, energy_exchange, plane_displacement)
                        h_bbbb_media += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide, i + 1, 0, j, 0, penalty_1, self._pair_energies, energy_phase_term))
                    except (IndexError, KeyError): pass
                
                    try:
                        energy_phase_term = Auxiliar_mj_modification.calculate_pairs_mj_modification(i, j-1, self._pair_energies_phase_1, self._pair_energies_phase_2, sign_function, weight_interface, energy_exchange, plane_displacement)
                        h_bbbb_media += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide, i, 0, j - 1, 0, penalty_1, self._pair_energies, energy_phase_term))
                    except (IndexError, KeyError): pass
                
                    try:
                        energy_phase_term = Auxiliar_mj_modification.calculate_pairs_mj_modification(i, j+1, self._pair_energies_phase_1, self._pair_energies_phase_2, sign_function, weight_interface, energy_exchange, plane_displacement)
                        h_bbbb_media += (self._contact_map.lower_main_upper_main[i][j]) ^ (self._distance_map.second_neighbor(self._peptide, i, 0, j + 1, 0, penalty_1, self._pair_energies, energy_phase_term))
                    except (IndexError, KeyError): pass
                
                    # Correct qubits if necessary
                    h_bbbb_media = _fix_qubits(h_bbbb_media, self._has_side_chain_second_bead)

        else:
            h_bbbb_media = 0

   
        return h_bbbb_media




    




    def _create_h_bbsc_and_h_scbb(self) -> Union[PauliSumOp, PauliOp]:
        """
        Creates Hamiltonian term corresponding to 1st neighbor interaction between main/backbone (BB) and side chain (SC) beads. 
        In the absence of side chains, this function returns a value of 0.

        Returns:
            Tuple of Hamiltonian terms consisting of backbone and side chain interactions.
        """
        penalty_1 = self._penalty_parameters.penalty_1
        h_bbsc = 0
        h_scbb = 0
        main_chain_len = len(self._peptide.get_main_chain)
        side_chain = self._peptide.get_side_chain_hot_vector()
        for i in range(1, main_chain_len - 3):
            for j in range(i + 4, main_chain_len + 1):
                if (j - i) % 2 == 1: continue

                if side_chain[j - 1] == 1:
                    h_bbsc += self._contact_map.lower_main_upper_side[i][j] ^ (
                        self._distance_map.first_neighbor(self._peptide, i, 0, j, 1, penalty_1, self._pair_energies)
                        + self._distance_map.second_neighbor(self._peptide, i, 0, j, 0, penalty_1, self._pair_energies))
                    try: h_bbsc += self._contact_map.lower_side_upper_side[i][j] ^ self._distance_map.first_neighbor(self._peptide, i, 1, j, 1, penalty_1, self._pair_energies)
                    except (IndexError, KeyError, TypeError): pass
                    try: h_bbsc += self._contact_map.lower_main_upper_side[i][j] ^ self._distance_map.second_neighbor(self._peptide, i + 1, 0, j, 1, penalty_1, self._pair_energies,)
                    except (IndexError, KeyError, TypeError): pass
                    try: h_bbsc += self._contact_map.lower_main_upper_side[i][j] ^ self._distance_map.second_neighbor(self._peptide,i - 1, 0, j, 1, penalty_1, self._pair_energies,)
                    except (IndexError, KeyError, TypeError): pass
                    
                if side_chain[i - 1] == 1:
                    h_scbb += self._contact_map.lower_side_upper_main[i][j] ^ (
                        self._distance_map.first_neighbor(self._peptide, i, 1, j, 0, penalty_1, self._pair_energies)
                        + self._distance_map.second_neighbor(self._peptide, i, 0, j, 0, penalty_1, self._pair_energies))
                    try: h_scbb += self._contact_map.lower_side_upper_main[i][j] ^ self._distance_map.second_neighbor(self._peptide, i, 1, j, 1, penalty_1, self._pair_energies)
                    except (IndexError, KeyError, TypeError): pass
                    try: h_scbb += self._contact_map.lower_side_upper_main[i][j] ^ self._distance_map.second_neighbor(self._peptide, i, 1, j + 1, 0, penalty_1, self._pair_energies,)
                    except (IndexError, KeyError, TypeError): pass
                    try: h_scbb += self._contact_map.lower_side_upper_main[i][j] ^ self._distance_map.second_neighbor(self._peptide, i, 1, j - 1, 0, penalty_1, self._pair_energies,)
                    except (IndexError, KeyError, TypeError): pass

        h_bbsc = _fix_qubits(h_bbsc, self._has_side_chain_second_bead)
        h_scbb = _fix_qubits(h_scbb, self._has_side_chain_second_bead)
        return h_bbsc, h_scbb

    def _create_h_scsc(self) -> Union[PauliSumOp, PauliOp]:
        """
        Creates Hamiltonian term corresponding to 1st neighbor interaction between side chain (SC) beads. 
        In the absence of side chains, this function returns a value of 0.

        Returns:
            Hamiltonian term consisting of side chain pairwise interactions
        """
        penalty_1 = self._penalty_parameters.penalty_1
        h_scsc = 0
        main_chain_len = len(self._peptide.get_main_chain)
        side_chain = self._peptide.get_side_chain_hot_vector()
        for i in range(1, main_chain_len - 3):
            for j in range(i + 5, main_chain_len + 1):
                if (j - i) % 2 == 0: continue
                if side_chain[i - 1] == 0 or side_chain[j - 1] == 0: continue
                h_scsc += self._contact_map.lower_side_upper_side[i][j] ^ (
                    self._distance_map.first_neighbor(self._peptide, i, 1, j, 1, penalty_1, self._pair_energies)
                    + self._distance_map.second_neighbor(self._peptide, i, 1, j, 0, penalty_1, self._pair_energies)
                    + self._distance_map.second_neighbor(self._peptide, i, 0, j, 1, penalty_1, self._pair_energies))
        return _fix_qubits(h_scsc, self._has_side_chain_second_bead)

    def _create_h_short(self) -> Union[PauliSumOp, PauliOp]:
        """
        Creates Hamiltonian constituting interactions between beads that are no more than 4 beads apart. If no side chains are present, this function returns 0.

        Returns:
            Contribution to energetic Hamiltonian for interactions between beads that are no more than 4 beads apart.
        """
        main_chain_len = len(self._peptide.get_main_chain)
        side_chain = self._peptide.get_side_chain_hot_vector()
        h_short = 0
        for i in range(1, main_chain_len - 2):
            # checks interactions between beads no more than 4 beads apart
            if side_chain[i - 1] == 1 and side_chain[i + 2] == 1:
                op1 = self._create_turn_operators(self._peptide.get_main_chain[i + 1], self._peptide.get_main_chain[i - 1].side_chain[0],)
                op2 = self._create_turn_operators(self._peptide.get_main_chain[i - 1], self._peptide.get_main_chain[i + 2].side_chain[0],)
                coeff = float(self._pair_energies[i][1][i + 3][1] + 0.1 * (self._pair_energies[i][1][i + 3][0] + self._pair_energies[i][0][i + 3][1]))
                composed = op1 @ op2
                h_short += (coeff * composed).reduce()
        h_short = _fix_qubits(h_short, self._has_side_chain_second_bead)

        return h_short
