# Quantum Computing in Peptide Folding Targeting Cellular Membranes

Modification of the IBM public code by adding protein folding in different phrases as well as the study of folding at a polar/apolar interface.


<p align="center">
  <img src="https://github.com/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes/assets/73825210/e2779b12-2022-427c-8b87-d5c72c4669c4" />
</p>


<div align="center">
  
[![GitHub forks](https://img.shields.io/github/forks/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)](https://github.com/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes/network)
[![GitHub issues](https://img.shields.io/github/issues/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)](https://github.com/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)
![GitHub pull requests](https://img.shields.io/github/issues-pr/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)
[![GitHub stars](https://img.shields.io/github/stars/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)](https://github.com/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes/stargazers)
![GitHub repo size](https://img.shields.io/github/repo-size/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)
[![GitHub license](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![GitHub language count](https://img.shields.io/github/languages/count/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes)
![Static Badge](https://img.shields.io/badge/Version-1.0-blue)

</div>

<div align="center">
  
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
[![made-with-bash](https://img.shields.io/badge/-Made%20with%20Bash-1f425f.svg?logo=image%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyZpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw%2FeHBhY2tldCBiZWdpbj0i77u%2FIiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8%2BIDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuNi1jMTExIDc5LjE1ODMyNSwgMjAxNS8wOS8xMC0wMToxMDoyMCAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENDIDIwMTUgKFdpbmRvd3MpIiB4bXBNTTpJbnN0YW5jZUlEPSJ4bXAuaWlkOkE3MDg2QTAyQUZCMzExRTVBMkQxRDMzMkJDMUQ4RDk3IiB4bXBNTTpEb2N1bWVudElEPSJ4bXAuZGlkOkE3MDg2QTAzQUZCMzExRTVBMkQxRDMzMkJDMUQ4RDk3Ij4gPHhtcE1NOkRlcml2ZWRGcm9tIHN0UmVmOmluc3RhbmNlSUQ9InhtcC5paWQ6QTcwODZBMDBBRkIzMTFFNUEyRDFEMzMyQkMxRDhEOTciIHN0UmVmOmRvY3VtZW50SUQ9InhtcC5kaWQ6QTcwODZBMDFBRkIzMTFFNUEyRDFEMzMyQkMxRDhEOTciLz4gPC9yZGY6RGVzY3JpcHRpb24%2BIDwvcmRmOlJERj4gPC94OnhtcG1ldGE%2BIDw%2FeHBhY2tldCBlbmQ9InIiPz6lm45hAAADkklEQVR42qyVa0yTVxzGn7d9Wy03MS2ii8s%2BeokYNQSVhCzOjXZOFNF4jx%2BMRmPUMEUEqVG36jo2thizLSQSMd4N8ZoQ8RKjJtooaCpK6ZoCtRXKpRempbTv5ey83bhkAUphz8fznvP8znn%2B%2F3NeEEJgNBoRRSmz0ub%2FfuxEacBg%2FDmYtiCjgo5NG2mBXq%2BH5I1ogMRk9Zbd%2BQU2e1ML6VPLOyf5tvBQ8yT1lG10imxsABm7SLs898GTpyYynEzP60hO3trHDKvMigUwdeaceacqzp7nOI4n0SSIIjl36ao4Z356OV07fSQAk6xJ3XGg%2BLCr1d1OYlVHp4eUHPnerU79ZA%2F1kuv1JQMAg%2BE4O2P23EumF3VkvHprsZKMzKwbRUXFEyTvSIEmTVbrysp%2BWr8wfQHGK6WChVa3bKUmdWou%2BjpArdGkzZ41c1zG%2Fu5uGH4swzd561F%2BuhIT4%2BLnSuPsv9%2BJKIpjNr9dXYOyk7%2FBZrcjIT4eCnoKgedJP4BEqhG77E3NKP31FO7cfQA5K0dSYuLgz2TwCWJSOBzG6crzKK%2BohNfni%2Bx6OMUMMNe%2Fgf7ocbw0v0acKg6J8Ql0q%2BT%2FAXR5PNi5dz9c71upuQqCKFAD%2BYhrZLEAmpodaHO3Qy6TI3NhBpbrshGtOWKOSMYwYGQM8nJzoFJNxP2HjyIQho4PewK6hBktoDcUwtIln4PjOWzflQ%2Be5yl0yCCYgYikTclGlxadio%2BBQCSiW1UXoVGrKYwH4RgMrjU1HAB4vR6LzWYfFUCKxfS8Ftk5qxHoCUQAUkRJaSEokkV6Y%2F%2BJUOC4hn6A39NVXVBYeNP8piH6HeA4fPbpdBQV5KOx0QaL1YppX3Jgk0TwH2Vg6S3u%2BdB91%2B%2FpuNYPYFl5uP5V7ZqvsrX7jxqMXR6ff3gCQSTzFI0a1TX3wIs8ul%2Bq4HuWAAiM39vhOuR1O1fQ2gT%2F26Z8Z5vrl2OHi9OXZn995nLV9aFfS6UC9JeJPfuK0NBohWpCHMSAAsFe74WWP%2BvT25wtP9Bpob6uGqqyDnOtaeumjRu%2ByFu36VntK%2FPA5umTJeUtPWZSU9BCgud661odVp3DZtkc7AnYR33RRC708PrVi1larW7XwZIjLnd7R6SgSqWSNjU1B3F72pz5TZbXmX5vV81Yb7Lg7XT%2FUXriu8XLVqw6c6XqWnBKiiYU%2BMt3wWF7u7i91XlSEITwSAZ%2FCzAAHsJVbwXYFFEAAAAASUVORK5CYII%3D)](https://www.gnu.org/software/bash/)

</div>



# Description of the challenge ⚛️:
<p align="justify"> 
Antimicrobial peptides (AMPs) are critical components of the innate immune system that are present in all living organisms. These peptides have been primarily associated with a defensive role against exogenous infections caused by bacteria, viruses and fungi, and are considered powerful and versatile endogenous antibiotics, capable of resisting bacterial adaptation for millions of years. However, recent research advances have pointed to the link between AMPs and a broader spectrum of diseases, such as cancer and various human inflammatory and autoimmune diseases, including ageing.
In general, the study of peptide and protein folding is an intrinsically very complex problem whose practical solution is beyond the reach of classical algorithms even in very simplified models such as the one that considers amino acids as simple spheres with different hydrophobic-polar character. In this scenario, quantum computers emerge as a promising tool despite the current state of intermediate noise level technology (NISQ). Recent work has attempted to solve this problem for relatively short amino acid sequences  in homogeneous media. In relation to our study, the complexity is further increased because the structure and functionality of AMPs are strongly influenced by the inhomogeneous and anisotropic environment in which they exhibit their activity, specifically, the membrane-water interface. This particularity highlights the need to intensify efforts in the development of quantum computing techniques, which may lead to significant discoveries in the study and development of new antimicrobial peptides, thus expanding the frontiers of our scientific understanding and therapeutic potential.
</p>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#introduction-movie_camera">Introduction</a></li>
    <li><a href="#general-pre-requirements-computer">Pre-requirements</a></li>
    <li><a href="#required-files-">Required-files</a></li>
    <li><a href="#usage-%EF%B8%8F">Usage</a></li>
    <li><a href="#output-files-">Output-files</a></li>
    <li><a href="#wiki-">Wiki</a></li>
    <li><a href="#distribution-of-tasks-%EF%B8%8F">Distribution of tasks</a></li>
    <li><a href="#faqs-interrobang">FAQs</a></li>
    <li><a href="#contributing-%EF%B8%8F">Contributing</a></li>
    <li><a href="#cite-mortar_board">Cite</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
    <li><a href="#license-green_book">License</a></li>
  </ol>
</details>

<!-- Introduction --> 
## Introduction :movie_camera:

<p align="justify"> 
In this repository, quantum computational methods have been developed and optimised to efficiently and accurately address the folding of amino acid sequences (peptides) at a polar/apolar interface, emulating a realistic model of the water/membrane interface, which is fundamental to understand the structure and function of Antimicrobial Peptides (AMPs) in specific biological contexts. It also allows the study of folding in both polar and non-polar monophases.
</p>

<!-- Pre-requirements -->
## General Pre-requirements :computer:

In order to run the program it is necessary to have the following requirements installed:

The program is written in Python language so Python version 3 or higher is required. Also, for python programs to work properly the following libraries are needed:
- [Qiskit](https://www.ibm.com/quantum/qiskit)
- [Matplotlib](https://matplotlib.org)
- [Numpy](https://numpy.org)
- [Re](https://docs.python.org/3/library/re.html)
- [Os](https://docs.python.org/3/library/os.html)
- [Mayavi](http://docs.enthought.com/mayavi/mayavi/)
- [Imageio](https://pypi.org/project/imageio/) 
- [Sys](https://docs.python.org/3/library/sys.html)
- [Time](https://docs.python.org/3/library/time.html)
- [Collections](https://docs.python.org/3/library/collections.html)
- [Glob](https://docs.python.org/3/library/glob.html)
- [Fileinput](https://docs.python.org/es/3/library/fileinput.html)


To be able to use these libraries, it is recommended to work with anaconda environment.


If you are working in an anaconda environment and you need to install a specific module:

```
conda install module
```

or, if you have pip installed:

```
conda install pip
```
```
pip install module
```


<!-- Required-files -->
## Required files 📋
No external file is required.

<!-- Usage -->
## Usage ⚙️
In order to run this program, the following command has to be used:
```
make run # To perform the search for the most stable state of your amino acid sequence
make representation # To get a picture of the most stable structure of your system.
make plots # It allows to obtain the evolution of the energy with the number of iterations.
```
### In the main.py file, in the inputs section we can choose:

> [!NOTE]  
  INPUTS:
  
    * main_chain: string            # amino acid sequence
    * ws_phase_1: float             # average interaction of an amino acid with solvent in phase 1
    * ws_phase_2: float             # average interaction of an amino acid with solvent in phase 2
    * cs_phase_1: float             # model solvent in phase 1
    * cs_phase_2: float             # model solvent in phase 2
    * exchange_solvent: float       # energy exchange between interfaces
    * correction_mj: bool           # indicates whether we decide to use the Miyazawa-Jerningan potential upgrade or not.
    * one_solvent_parameter: bool   # These options do not affect the calculation of two phases, is only if we want to add a parameter to modulate the interaction of amino acids according to their hydrophobicity value with the phase.
    * penalty_back: float           # A penalty parameter used to penalize turns along the same axis. This term is used to eliminate sequences where the same axis is chosen twice in a row. In this way we do not allow for a chain to fold back into itself.
    * penalty_chiral: float         # A penalty parameter used to impose the right chirality.
    * penalty_1: float              # A penalty parameter used to penalize local overlap between beads within a nearest neighbor contact.
    * maxiter: float                # Number of iterations to be performed.


Not interface when: 

```
cs_phase_1 == cs_phase_2 and ws_phase_1 == ws_phase_2
```

The original IBM work can be obtained by setting up:

```
* ws_1 = ws_2
* cs_1 = cs_2
* exchange_sovent_solvent = 0
* correction = False
* one_solvent_paramenter = False
* penalty_back = 10
* penalty_chiral = 10
* penalty_1 = 10
```


### In the qiskit_research/protein_folding/interface_bead_interaction/interface_parameters.py file:
> [!IMPORTANT]
Allows us to modify different parameters related to the creation of the interface.

``` 
* displacement_of_the_interface_plane: float  # Value for the displacement of the interface plane. 
* axis: int                                   # Axis perpendicular to the interface.  Axis to be selected for defining the plane of the interface: The possible values are 0, 1, 2 or 3.
* cls._instance._weight_interface: abs(float) # Value for the importance of the interface. Initial value must be positive!
```             

Value for the weight to be given to the interface Hamiltonian:

* If weight_interface > 0, the hydrophobic solvent goes below the plane and the hydrophilic solvent goes on top.
* If weight_interface < 0, the hydrophobic solvent  goes on top of the plane and the hydrophilic solvent goes below.

The sign will be given according to the selected values of phase 1 cs and phase 2 cs: 

```
if cs_phase_1 > cs_phase_2:
    mj_modified_parameters = -1
elif cs_phase_2 > cs_phase_1:
    mj_modified_parameters = 1
```

Therefore, in the first case the importance of the interface will be ``` cls._instance._weight_interface * -1 ``` and in the second case ``` cls._instance._weight_interface * 1 ```.


### Quantum simulation parameters
> [!TIP]
Brief summary:
* cs > 0: model hydrophilic solvents (water)
* cs < 0: model hydrophobic solvents (oil)
* w: describes the average interaction of an amino acid with solvent
* w > 0: aminoacids repel the solvent
* w < 0:  aminoacids interact on average attractively with the solvent


<!-- Output files -->
## Output files 📋
* **Data_plane.txt:**  Contains the equation of the plane of the interface of the form Ax+By+Cz+D=0.

* **energy_iteractions.txt:** Presents information on the evolution of energy (VQE method), with respect to each interaction.

* **output.xyz:** File that presents the conformations of the peptide that I explore during the calculation. Can be loaded with [VMD](https://www.ks.uiuc.edu/Research/vmd).

* **plane_vmd.xyz:** File that presents the interface in a .xyz file. Can be loaded with [VMD](https://www.ks.uiuc.edu/Research/vmd).

* **STABLE.xyz:** File that presents the most stable conformation of the peptide (minimum energy). Can be loaded with [VMD](https://www.ks.uiuc.edu/Research/vmd).

<p align="center">
  <img src="https://github.com/DanielCondeTorres/Quantum-Quirks-Unraveling-Protein-Folding-Through-an-Advanced-Interface./assets/117435891/456188df-bc99-4fa7-adcd-7449412f5ee3" />
</p>


* **Results.txt:**  File that stores information about the computing time and the parameters used in your simulation.



<!-- WIKI -->
## Wiki 📖

Additional information about this modifications can be found in our [Wiki](https://github.com/TeamMduse/Quantum-Computing-in-Peptide-Folding-Targeting-Cellular-Membranes/wiki) **(...Work in progress...)**.


<!-- DISTRIBUTION OF TASKS -->
## Distribution of tasks ✒️ 
Project coordinator: Ándres Gómez, Ángel Piñeiro and Rebeca García-Fandino

- Main program (Python program): Daniel Conde-Torres, Mariamo Mussa-Juane, Daniel Faílde, Ángel Piñeiro, Ándres Gómez, Ángel Piñeiro and Rebeca García-Fandino
- GitHub Designer: Daniel Conde-Torres, Mariamo Mussa-Juane and Daniel Faílde


<!-- FAQs -->
## FAQs :interrobang:


### If I don't have a quantum computer, can I use it?
Yes, it is designed to run on quantum computer emulators.

### What sequence size should I use on a conventional computer? 
In a reasonable time, results can be obtained with sequences of 8 amino acids (104.63396787643 s). However, with a sequence of 10 amino acids you would need 102.3 h!


<!-- CONTRIBUTING -->
## Contributing 🖇️
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

<!-- Cite --> 
## Cite :mortar_board:

<!-- Acknowledgements --> 
## Acknowledgements:

This work was supported by the Spanish Agencia Estatal de Investigación (AEI) and the ERDF (PID2022-141534OB-I00 and PID2019-111373RB-I00), by Xunta de Galicia and the ERDF (02\_IN606D\_2022\_2667887, ED431F 2020/05, ED431C 2021/21, ED431B 2022/36 and Centro singular de investigación de Galicia accreditation 2016-2019, ED431G/09, Ministerio de Universidadese for the predoctoral contract (FPU22/00636) and the REACT-EU Funds. All calculations were carried out at CESGA.

<!-- License --> 
## License :green_book:
Apache License
Version 2.0, January 2004
http://www.apache.org/licenses/
