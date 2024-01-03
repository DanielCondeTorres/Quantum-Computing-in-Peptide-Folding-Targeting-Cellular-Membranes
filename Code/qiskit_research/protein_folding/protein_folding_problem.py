# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Defines a protein folding problem that can be passed to algorithms."""
from __future__ import annotations

from typing import TYPE_CHECKING, List, Union

from qiskit.algorithms.minimum_eigensolvers import MinimumEigensolverResult
from qiskit.opflow import PauliOp, PauliSumOp

from .interactions.interaction import Interaction
from .penalty_parameters import PenaltyParameters
from .peptide.peptide import Peptide
from .qubit_op_builder import QubitOpBuilder
from .qubit_utils import qubit_number_reducer
from .sampling_problem import SamplingProblem

if TYPE_CHECKING:
    from .protein_folding_result import ProteinFoldingResult


class ProteinFoldingProblem(SamplingProblem):
    """Defines a protein folding problem that can be passed to algorithms. Example initialization:

    .. code-block:: python

        penalty_terms = PenaltyParameters(15, 15, 15)
        main_chain_residue_seq = "SAASSASAAG"
        side_chain_residue_sequences = ["", "", "A", "A", "A", "A", "A", "A", "S", ""]
        peptide = Peptide(main_chain_residue_seq, side_chain_residue_sequences)
        mj_interaction = MiyazawaJerniganInteraction()
        protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
        qubit_op = protein_folding_problem.qubit_op()
    """

    def __init__(
        self,
        peptide: Peptide,
        interaction: Interaction,
        penalty_parameters: PenaltyParameters,
        phase_w_1: float, 
        phase_cs_1: float, 
        phase_w_2:float,
        phase_cs_2:float,
        exchange_solvent_solvent:float,
        correction_mj: bool,
        one_solvent_parameter: bool
    ):
        """
        Args:
            * peptide: A peptide object that defines the protein subject to the folding problem.
            * interaction: A type of interaction between the beads of the peptide.
            * penalty_parameters: Parameters that define the strength of constraints enforcing in
                                the problem.
            * ws_phase_1: average interaction of an amino acid with solvent in phase 1
            * ws_phase_2: average interaction of an amino acid with solvent in phase 2
            * cs_phase_1: model solvent in phase 1
            * cs_phase_2: model solvent in phase 2
            * exchange_solvent: energy exchange between interfaces
            * correction_mj: indicates whether we decide to use the Miyazawa-Jerningan potential upgrade or not.
            * one_solvent_parameter: These options do not affect the calculation of two phases, is only if we want to add a parameter to modulate the interaction of amino acids according to their hydrophobicity value with the phase.   
        """
        self._peptide = peptide
        self._interaction = interaction
        self._penalty_parameters = penalty_parameters

        self._pair_energies = interaction.calculate_energy_matrix(peptide.get_main_chain.main_chain_residue_sequence,)
        

        self._phase_w_1 = phase_w_1
        self._phase_w_2 = phase_w_2
        self._phase_cs_1 = phase_cs_1
        self._phase_cs_2 = phase_cs_2
        self._one_solvent_parameter = one_solvent_parameter
        # We apply the correction to the Miyazawa-Jerningan potential if this option has been chosen.
        self._pair_energies_correction_mj = interaction.calculate_media_energy_matrix(peptide.get_main_chain.main_chain_residue_sequence,
                                                                                  phase_w_1, 
                                                                                  phase_cs_1, 
                                                                                  phase_w_2,
                                                                                  phase_cs_2,
                                                                                  exchange_solvent_solvent,
                                                                                  correction_mj
                                                                                     )
        self._correction_mj = correction_mj
        self._qubit_op_builder = QubitOpBuilder(self._peptide, self._pair_energies, self._penalty_parameters,self._pair_energies_correction_mj, self._phase_cs_1, self._phase_cs_2,self._phase_w_1, self._phase_w_2, self._one_solvent_parameter,self._correction_mj )

        self._unused_qubits: List[int] = []
    

    def qubit_op(self) -> Union[PauliSumOp, PauliOp]:
        """
        Builds a qubit operator for the Hamiltonian encoding a protein folding problem. The
        number of qubits needed for optimization is optimized (compressed), if possible.
        To obtain the full qubit operator for a Hamiltonian, use the method `qubit_op_full`.

        Returns:
            A qubit operator for the Hamiltonian encoding a protein folding problem on an
            optimized number of qubits.
        """
        qubit_operator, unused_qubits = qubit_number_reducer.remove_unused_qubits(
            self._qubit_op_full()
        )
        self._unused_qubits = unused_qubits
        return qubit_operator

    def _qubit_op_full(self) -> Union[PauliOp, PauliSumOp]:
        """
        Builds a full qubit operator for the Hamiltonian encoding a protein folding problem. Full
        means that the number of qubits needed for optimization is not optimized and may be
        larger that necessary. To ensure the optimal number of qubits, use the method `qubit_op`.

        Returns:
            A qubit operator for the Hamiltonian encoding a protein folding problem.
        """
        qubit_operator = self._qubit_op_builder.build_qubit_op()
        return qubit_operator

    def interpret(self, raw_result: MinimumEigensolverResult) -> "ProteinFoldingResult":
        """
        Interprets the raw algorithm result, in the context of this problem, and returns a
        ProteinFoldingResult. The returned class can plot the protein and generate a
        .xyz file with the coordinates of each of its atoms.
        Args:
            raw_result: The raw result of solving the protein folding problem.

        Returns:
            A :class:`~qiskit_research.protein_folding.ProteinFoldingResult`
            instance that contains the protein folding result.
        """
        # pylint: disable=import-outside-toplevel
        from .protein_folding_result import ProteinFoldingResult

        probs = raw_result.eigenstate.binary_probabilities()
        best_turn_sequence = max(probs, key=probs.get)
        return ProteinFoldingResult(
            unused_qubits=self.unused_qubits,
            peptide=self.peptide,
            turn_sequence=best_turn_sequence,
        )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def interpret_intermediate(self, bitstring: MinimumEigensolverResult) -> "ProteinFoldingResult":
        """
        Interprets the raw algorithm intermediates, in the context of this problem, and returns a
        ProteinFoldingResult. The returned class can plot the protein and generate a
        .xyz file with the coordinates of each of its atoms.
        
        Args:
            *bitstring: The raw conformations

        Returns:
            A :class:`~qiskit_research.protein_folding.ProteinFoldingResult`
            instance that contains the protein folding conformation of the bitstring.
        """
        from .protein_folding_result import ProteinFoldingResult
        return ProteinFoldingResult(
            unused_qubits=self.unused_qubits,
            peptide=self.peptide,
            turn_sequence=bitstring,
        )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @property
    def unused_qubits(self) -> List[int]:
        """Returns the list of indices for qubits in the original problem formulation that were
        removed during compression."""
        return self._unused_qubits

    @property
    def peptide(self) -> Peptide:
        """Returns the peptide defining the protein subject to the folding problem."""
        return self._peptide
