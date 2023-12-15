#---------------------------------------#
#-------------Libraries-----------------#
#---------------------------------------#
# Now, using qiskit-research:
import sys
sys.path.append("..") # Add previous path to the sistem
import time
from qiskit_research.protein_folding.protein_folding_problem import (ProteinFoldingProblem)
from qiskit_research.protein_folding.interactions.random_interaction import (RandomInteraction)
from qiskit_research.protein_folding.interactions.miyazawa_jernigan_interaction import (MiyazawaJerniganInteraction)
from qiskit_research.protein_folding.peptide.peptide import Peptide
from qiskit_research.protein_folding.protein_folding_problem import (ProteinFoldingProblem)
from qiskit_research.protein_folding.penalty_parameters import PenaltyParameters
import matplotlib.pyplot as plt
from qiskit.utils import algorithm_globals, QuantumInstance
#Using VQE with CVaR expectation value for the solution of the problem
from qiskit.circuit.library import RealAmplitudes
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.algorithms.minimum_eigensolvers import SamplingVQE
from qiskit import execute, Aer
from qiskit.primitives import Sampler
import os
from collections import Counter
import fileinput
import glob
import re



class VQEHandler:
    ''' Modification of the previous IBM code function to get the bitstring in each iteration.'''
    def __init__(self,ansatz):
        # We initialise a series of empty lists to store results of interest.     
        self.counts = []
        self.values = []
        self.bit_strings = []
        self.q_inicial = []
        self.ansatz = ansatz
        self.energy_values = []
   
    def store_intermediate_result(self, eval_count, parameters, mean, std):
        self.counts.append(eval_count)
        self.values.append(mean)
        # Configurate ansatz parameters
        parameter_dict = dict(zip(self.ansatz.parameters, parameters))
        circuit = self.ansatz.bind_parameters(parameter_dict)

        # Quantum state measurement
        circuit.measure_all()

        # Run the circuit in a simulator
        backend = Aer.get_backend('qasm_simulator')
        job = execute(circuit, backend, shots=1)
        result = job.result()
        
        #  Get backend configuration
        backend_conf = backend.configuration()

        # Get the number of qubits in the configuration
        num_qubits = backend_conf.n_qubits

        self.q_inicial.append(num_qubits)
        print(f"Number of qubits available in the simulator: {num_qubits}")

        counts_p = result.get_counts()
        bitstring = list(counts_p.keys())[0]
        print(f"Iteration {eval_count}: energy = {mean}, bitstring = {bitstring}")
        self.bit_strings.append(bitstring)
        self.energy_values.append(mean)




#The Protein consists of a main chain that is a linear chain of aminoacids.
def search_for_the_most_stable_conformation(main_chain : str, ws_phase_1 : float, cs_phase_1: float, ws_phase_2 : float, cs_phase_2:  float,  exchange_solvent_solvent: float = 0, correction_mj: bool = True, reps: int = 50, maxiter: int = 1, one_solvent_parameter: bool = False,  penalty_back = 10,  penalty_chiral = 10, penalty_1 = 10):
    '''
    
    
    This function generates different graphs and files for a graphical representation:

    INPUTS:
   
    * main_chain: amino acid sequence
    * ws_phase_1: average interaction of an amino acid with solvent in phase 1
    * ws_phase_2: average interaction of an amino acid with solvent in phase 2
    * cs_phase_1: model solvent in phase 1
    * cs_phase_2: model solvent in phase 2
    * exchange_solvent: energy exchange between interfaces
    * correction_mj: indicates whether we decide to use the Miyazawa-Jerningan potential upgrade or not.
    * reps: 
    * maxiter: 
    * one_solvent_parameter: These options do not affect the calculation of two phases, is only if we want to add a parameter to modulate the interaction of amino acids according to their hydrophobicity value with the phase.
    * penalty_back: A penalty parameter used to penalize turns along the same axis. This term is used to eliminate sequences where the same axis is chosen twice in a row. In this way we do not allow for a chain to fold back into itself.
    * penalty_chiral: A penalty parameter used to impose the right chirality.
    * penalty_1: A penalty parameter used to penalize local overlap between beads within a nearest neighbor contact.




    OUTPUTS:


    '''
    #Beyond the main chain of the protein there may be aminoacids attached to the residues of the main chain. Our model allows for side chains of the maximum length of one. Elongated side chains would require the introduction of additional penalty terms which are still under development. In this example we do not consider any side chains to keep the real structure of the neuropeptide.
    side_chains = [""] * len(main_chain)
    
    #Interacion between Aminoacids
    random_interaction = RandomInteraction()#Beyond this model we also allow for random contact maps (interactions) that provide a random interaction map. One can also introduce a custom interaction map that enhances certain configurations of the protein (e.g. alpha helix, beta sheet etc).
    mj_interaction = MiyazawaJerniganInteraction()
    
    #Physical Constraints to ensure that physical constraints are respected
    penalty_terms = PenaltyParameters(penalty_chiral, penalty_back, penalty_1)
    
    #Peptide definition
    peptide = Peptide(main_chain, side_chains)#Based on the main chain and possible side chains we define the peptide object that includes all the structural information of the modeled system.
    
    #Protein folding problem
    protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms,ws_phase_1,cs_phase_1,ws_phase_2,cs_phase_2,exchange_solvent_solvent, correction_mj, one_solvent_parameter)#Based on the defined peptide, the interaction (contact map) and the penalty terms we defined for our model we define the protein folding problem that returns qubit operators.
    qubit_op = protein_folding_problem.qubit_op()
    
    # set classical optimizer:  COBYLA for the classical optimization part
    optimizer = COBYLA(maxiter=maxiter)
    # set variational ansatz
    ansatz = RealAmplitudes(reps=reps)    
    
    #  Quantum variational algorithm
    vqe_handler = VQEHandler(ansatz) 
    vqe = SamplingVQE(Sampler(),
                      ansatz=ansatz,
                      optimizer=optimizer,
                      aggregation=0.1,
                      callback=vqe_handler.store_intermediate_result)
    
    # Result of the problem
    raw_result = vqe.compute_minimum_eigenvalue(qubit_op)
    result = protein_folding_problem.interpret(raw_result=raw_result)
    
    return result, vqe_handler.counts, vqe_handler.values, ansatz, peptide, vqe_handler.q_inicial, vqe_handler.bit_strings, protein_folding_problem, vqe_handler.energy_values



def extract_number(filename):
    
    '''
    This function is designed to extract the numeric portion from a filename. Here's the breakdown:
    '''
    
    # This line uses Python's regular expression module, re, to search for one or more digits (\d+) in the filename. The match object will store the first sequence of digits found in the filename.
    match = re.search(r'\d+', filename)

    # This line checks if any numeric sequence was found (if match). If found, it extracts the matched string using the group() method and then converts it to an integer using int(). If no match is found, it returns 0.
    return int(match.group()) if match else 0



def create_figure(protein_folding_problem, bitstring: list, counts: list, values: list, result, ansatz, main_chain: str, ws_phase_1: float, cs_phase_1: float, ws_phase_2: float, cs_phase_2: float, exchange_solvent: float, correction_mj: bool, one_solvent_parameter: bool,video: bool = False):

    '''
    This function generates different graphs and files for a graphical representation:  

    INPUTS:

    * protein_folding_problem:
    * bitstring: list of all conformations in bitstrng format
    * counts: 
    * values: 
    * result: 
    * ansatz: 
    * main_chain: amino acid sequence
    * ws_phase_1: average interaction of an amino acid with solvent in phase 1
    * ws_phase_2: average interaction of an amino acid with solvent in phase 2
    * cs_phase_1: model solvent in phase 1
    * cs_phase_2: model solvent in phase 2
    * exchange_solvent: energy exchange between interfaces
    * correction_mj: indicates whether we decide to use the Miyazawa-Jerningan potential upgrade or not.
    * video: if we want to save the .xyz video
    
    Most of these inputs are just to name the files and identify them later.


    OUTPUTS: 

    * output.xyz: generates an .xyz file in video format of all scanned conformations.
    * STABLE.xyz: shows the most stable conformation of our peptide.
    * Data_plane.txt: will show the equation of the plane of the interface, this file is not explicitly generated in this function.
    * Energy_ *.png: shows the variation of the energy with the iterations performed.
    * Conformation_*.png: shows the most stable conformation of our peptide without representing the plane. 

    '''
    
    #Valores organizar proyecto:
    plt.tight_layout()
    displacement_of_the_interface_plane = 0.5
    h_interface = 0.1
    eje = 1

    directorio =  'NEW'+str(main_chain)+'_'+'media1_'+str(cs_phase_1)+'media2_'+str(cs_phase_2)+'displacement'+str(displacement_of_the_interface_plane)+'h_interfece'+str(h_interface)+'eje'+str(eje)+str(correction_mj)+str(one_solvent_parameter)
     # This segment checks if a directory exists (based on the main_chain parameter) and creates one if it doesn't. Then, it saves the earlier generated plot with a specific naming convention. 
    if not os.path.exists(directorio):
        os.makedirs(directorio)

    # Plot Energy vs. Iterations:
    fig = plt.figure()
    plt.plot(counts, values)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Conformation Energy",fontsize=20)
    plt.xlabel("VQE Iterations",fontsize=20)
    fig.add_axes([0.44, 0.51, 0.44, 0.32])
    plt.plot(counts[40:], values[40:])
    with open(directorio+'/'+'energy_iteractions.txt', 'w') as file:
        for row in range(len(counts)):
            file.write(f"{counts[row]}\t\t{values[row].real}\n")
    plt.ylabel("Conformation Energy", fontsize=18)
    plt.xlabel("VQE Iterations",fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    # Save Energy Plot:
    
    plt.savefig(directorio+'/'+'Energy_'+str(main_chain)+'_'+'media1_'+str(ws_phase_1)+'_'+str(cs_phase_1)+'media2_'+str(ws_phase_2)+'_'+str(cs_phase_2)+'exchange_'+str(exchange_solvent)+str(correction_mj)+'.png')#, bbox_inches='tight')
    
    # Print Bitstring Representation:
    print("The bitstring representing the shape of the protein during optimization is: ",result.turn_sequence,)
    print(f"The folded protein's main sequence of turns is: {result.protein_shape_decoder.main_turns}")
    print(f"and the side turn sequences are: {result.protein_shape_decoder.side_turns}")
    
    # Generate and Save Protein Structure Plot:
    fig = result.get_figure(title="Protein Structure", ticks=False, grid=True)
    fig.get_axes()[0].view_init(10, 70)
    plt.savefig(directorio+'/'+'Conformation_'+str(main_chain)+'_'+'media1_'+str(ws_phase_1)+'_'+str(cs_phase_1)+'media2_'+str(ws_phase_2)+'_'+str(cs_phase_2)+'exchange_'+str(exchange_solvent)+str(correction_mj)+'.png')
    
    # Save Most Stable Conformation: This section tries to save the most stable conformation in an .xyz file. If a file named STABLE already exists, it deletes that file and saves anew.
    try:
        result.save_xyz_file(name='STABLE', path=directorio)
    except FileExistsError:
        os.remove(os.path.join(directorio, 'STABLE.xyz'))
        result.save_xyz_file(name='STABLE', path=directorio)

    # Generate and Save Video: When the video parameter is set to True, this block saves each protein conformation as individual .xyz files. 
    # It then combines these files into a single output.xyz file, creating a continuous representation of the protein's conformation changes. Any individual .xyz files apart from output.xyz and STABLE.xyz are then removed.
    if video == True:
        try:
            for elemento in range(len(bitstring)):
                result_intermediate = protein_folding_problem.interpret_intermediate(bitstring = bitstring[elemento])
                result_intermediate.save_xyz_file(name=str(elemento),path=directorio)
            #directorio = str(main_chain)+'_'+'media1_'+str(ws_phase_1)+'_'+str(cs_phase_1)+'media2_'+str(ws_phase_2)+'_'+str(cs_phase_2)+'exchange_'+str(exchange_solvent)+str(correction_mj)
            filenames = sorted(glob.glob(os.path.join(directorio, "*.xyz")), key=extract_number)
            with open(os.path.join(directorio, "output.xyz"), "w") as outfile:
                for filename in filenames:
                    with open(filename, "r") as infile:
                        for line in infile:
                            outfile.write(line)
            archivos_xyz = [archivo for archivo in os.listdir(directorio) if archivo.endswith(".xyz")]
            for archivo in archivos_xyz:
                if archivo not in [ "output.xyz","STABLE.xyz","plane_vmd.xyz"]:
                    archivo_path = os.path.join(directorio, archivo)
                    os.remove(archivo_path)
        except FileExistsError:
            pass
    else:
        pass

    # Return Statement:
    path = directorio
    return result.turn_sequence,path



from collections import Counter
import matplotlib.pyplot as plt

def histogram(bitstring: list, energies: list, path: str, top_n: float = 5, top_bars: int = 50):
    '''
    Generate a histogram with the frequency at which each bitstring appears, 
    to the top_n highest frequency values we add their energy value.
    
    * bitstring: list of all conformations in bitstring format
    * energies: list of energy values for each bitstring
    * path: place where we are going to save the file
    * top_n: how many of the higher frequency values we are going to consider to represent energy
    * top_bars: number of bars with the highest frequencies to represent on the histogram
    '''
    with open(path+'/'+'energy_frequenciues.txt', 'w') as file:
        for row in range(len(energies)):
            file.write(f"{energies[row].real}\t\t{bitstring[row]}\n")


    # Process the Bitstring List:
    result = [item.split()[0] for item in bitstring]

    # Calculate the Frequencies of each value in "result" and Percentages:
    total_items = len(result)
    frequencies = Counter(result)
    frequencies_percentage = {key: (value / total_items) * 100 for key, value in frequencies.items()}

    # Sort the frequencies and take the top 50:
    top_50_frequencies = dict(sorted(frequencies_percentage.items(), key=lambda item: item[1], reverse=True)[:top_bars])

    # Associate Bitstrings with Energies:
    bitstring_to_energy_new = dict(zip(result, energies))
    
    # Sort Frequencies and Determine Top Bitstrings:
    sorted_frequencies = sorted(top_50_frequencies.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_bitstrings = [item[0] for item in sorted_frequencies]

    # Plot the Histogram:
    plt.figure(figsize=(12, 7))
    bars = plt.bar(top_50_frequencies.keys(), top_50_frequencies.values(), color='lightcoral')

    # Label Top Bitstrings with Energy Values:
    for bitstring, freq in top_50_frequencies.items():
        if bitstring in top_bitstrings:
            energy = bitstring_to_energy_new[bitstring]
            height = freq
            plt.text(bitstring,
                     height,  # Adjusting the position for better visibility
                     f"{energy.real:.2f}",
                     ha='center',
                     color='blue',
                     fontsize=12,
                     rotation=45)

    # Additional Plot Customizations:
    plt.xlabel("Bitstrings",fontsize = 20)
    plt.ylabel("Frequency (%)", fontsize = 20)
    plt.xticks(rotation=90,fontsize=18)
    plt.yticks(fontsize = 18)
    plt.grid(axis='y')
    plt.tight_layout()

    # Save and Show Plot:
    plt.savefig(path + '/' + 'Modified_Histogram')
    plt.show()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                    Main program                                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
start_time = time.time() # Safe the time
algorithm_globals.random_seed = 23


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#            INPUTS          #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
main_chain = "WRDWGSGWDR" # Peptide sequence

'''

Modification of the Miyazawa- Jerningan potential following: 
Three-Dimensional Lattice Monte Carlo Simulations of Model Proteins. IV. Proteins at an Oil−Water Interface 
K. Leonhard, J. M. Prausnitz, and C. J. Radke Langmuir 2006 22 (7), 3265-3272 
DOI: 10.1021/la052535h 

Brief summary:

* cs > 0: model hydrophilic solvents (water)
* cs < 0: model hydrophobic solvents (oil)

* w: describes the average interaction of an amino acid with solvent
* w > 0: aminoacids repel the solvent
* w < 0:  aminoacids interact on average attractively with the solvent

* Exchange energy: is the phase 1 - phase 2 interaction

At an oil-water interface, typical values for these parameters are:

cs_water = 0.1 - 0.4
w_water = 0.16

cs_oil = -0.1
w_oil = In the study they change this parameter in the range: -0.04 - 0.36

Eoil_water = 0.3 - 0.5

'''

'''
Not interface when: cs_phase_1 == cs_phase_2 and ws_phase_1 == ws_phase_2
The original IBM work can be obtained by setting up:
* ws_1 = ws_2
* cs_1 = cs_2
* exchange_sovent_solvent = 0
* correction = False
* one_solvent_paramenter = False
'''
# Values of the first fase
ws_1 = 0.2
cs_1 = 0.8

# Values of the second fase
ws_2 = 0.2 
cs_2 = -0.4

# Energy exchange between interfaces
exchange_solvent_solvent = 0.

#This parameter allows you to indicate whether you want to apply the Miyazawa - Jerningan potential correction ( True ), or not (False). 
correction = False
# One_solvent_parameter: These options do not affect the calculation of two phases, is only if we want to add a parameter to modulate the interaction of amino acids according to their hydrophobicity value with the phase.
one_solvent_parameter = False
# Quantum simulation parameters
reps = 1 
maxiter = 200

'''
Note that in the folder: interface_bead_interaction/interface_parameters.py there are more (more subtle) parameters with which the user can also modify his system.
'''

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#      END   INPUTS          #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!#


# We perform quantum simulation
result,counts,values,ansatz,peptide,qubits_availables,bitstring,protein_folding_problem,energy = search_for_the_most_stable_conformation(main_chain, ws_1,cs_1, ws_2, cs_2, exchange_solvent_solvent, correction, reps, maxiter, one_solvent_parameter)

# We save the time it took to run the code
end_time = time.time() 
execution_time = end_time - start_time 



# Save information
information_file='Results.txt'
if not os.path.exists(information_file):
    # Create the header
    header = "#{:<20}\t\t{:<10}\t\t{:<5}\t\t{:<5}\t\t{:<5}\t\t{:<5}  \t{:<1}\t{:<5}\t{:<5}\t{:<5}  \t{:<1}\n".format(
                  "TIME(s)", "SEQUENCE", "WS_1", "CS_1", "WS_2", "CS_2", "EX_SOLVENT", "Reps_ansazt", "Iteracciones","Qb_availables", "Qb_used")
    with open(information_file, 'w') as information:
        information.write(header)
else:
        print(f" The file {information_file} already exists existe. No action taken.")
# Safe information in Results.txt
with open(information_file, 'a') as file: 
    file.write(f"{execution_time:<20.11f}\t\t{str(main_chain):<10}\t\t{ws_1:<5}\t\t{cs_1:<5}\t\t{ws_2:<5}\t\t{cs_2:<5}\t\t{exchange_solvent_solvent:<5}\t\t{reps:<5}\t\t{maxiter:<5}\t\t{qubits_availables[-1]:<5}\t\t{len(str(bitstring))}\t\tmod\n")

# Create a figure of the peptide
figure,path = create_figure(protein_folding_problem, bitstring,counts,values,result,ansatz,main_chain,ws_1,cs_1,ws_2,cs_2,exchange_solvent_solvent, correction,one_solvent_parameter, True)

# Create a histogram of all the bitstrings searched for.
histogram = histogram(bitstring[:],energy[:],path)

print(f"The execution time: {execution_time} s")
print(" Path : ", path)
