# In this class we include the values to touch for the calculation of the Hamiltonian of the interface, in the QubitOp class.

class interface_parameters:

    ''' In this class, the parameters to be used in the calculation of the Hamiltonian of the interface are defined, being free to be modified by the user. '''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ''' We created a dictionary of hydrophobicity values by following:
    Fauchère, J., and Pliska, V. 1983. Hydrophobic parameters {pi} of amino-acid side chains from the partitioning of N-acetyl-amino-acid amides. Eur. J. Med. Chem. 8: 369–375: 
    https://www.researchgate.net/publication/246404378_Hydrophobic_parameters_II_of_amino_acid_side-chains_from_the_partitioning_of_N-acetyl-amino_acid_amides '''

    AMINO = {'ASP': {'label':'D','hydrophobicity': -0.77,'behavior': 'acidic'},
             'GLU': {'label':'E','hydrophobicity': -0.64,'behavior': 'acidic'},
             'LYS': {'label':'K','hydrophobicity': -0.99,'behavior': 'basic'},
             'ARG': {'label':'R','hydrophobicity': -1.01,'behavior': 'basic'},
             'HIS': {'label':'H','hydrophobicity': 0.13,'behavior': 'basic'},
             'HISH': {'label':'H','hydrophobicity': 0.13,'behavior': 'basic'},
             'GLY': {'label':'G','hydrophobicity': -0.0,'behavior': 'hydrophobic'},
             'ALA': {'label':'A','hydrophobicity': 0.31,'behavior': 'hydrophobic'},
             'VAL': {'label':'V','hydrophobicity': 1.22,'behavior': 'hydrophobic'},
             'LEU': {'label':'L','hydrophobicity': 1.70,'behavior': 'hydrophobic'},
             'ILE': {'label':'I','hydrophobicity': 1.80,'behavior': 'hydrophobic'},
             'PRO': {'label':'P','hydrophobicity': 0.72,'behavior': 'hydrophobic'},
             'PHE': {'label':'F','hydrophobicity': 1.79,'behavior': 'hydrophobic'},
             'MET': {'label':'M','hydrophobicity': 1.23,'behavior': 'hydrophobic'},
             'TRP': {'label':'W','hydrophobicity': 2.25,'behavior': 'hydrophobic'},
             'SER': {'label':'S','hydrophobicity': -0.04,'behavior': 'polar'},
             'THR': {'label':'T','hydrophobicity': 0.26,'behavior': 'polar'},
             'CYS': {'label':'C','hydrophobicity': 1.54,'behavior': 'polar'},
             'TYR': {'label':'Y','hydrophobicity': 0.96,'behavior': 'polar'},
             'ASN': {'label':'N','hydrophobicity': -0.60,'behavior': 'polar'},
             'GLN': {'label':'Q','hydrophobicity': -0.22,'behavior': 'polar'}}

    ''' This website: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
    It includes several tables, so it can be modified. '''

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Igual aproximacion polinomica y no de Chebyshev, estamos seguros 
    ''' Coefficients of the Chebyshev polynomial for the sign function. '''
    coefficients_of_chebysev={'coef_7' : -0.00000312*0.5,
                              'coef_6' : 0.0,
                              'coef_5' : 0.00059*0.5,
                              'coef_4' : 0.0,
                              'coef_3' : -0.0364*0.5,
                              'coef_2' : 0.0,
                              'coef_1' : 0.9635*0.5,
                              'coef_0' : 0.0}
    ''' Each coef_i, i indicates the degree of the variable it accompanies. '''



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ''' Axis to be selected for defining the plane of the interface: The possible values are 0, 1, 2 or 3. '''

    axis = 1 #1
    '''  Axis perpendicular to the interface. '''

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ''' Value for the displacement of the interface plane. '''
    displacement_of_the_interface_plane = 1. #falta hacer despla 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ''' Value for the weight to be given to the interface Hamiltonian:
    
        * If weight_interface > 0, the hydrophobic solvent goes below the plane and the hydrophilic solvent goes on top.
        * If weight_interface < 0, the hydrophobic solvent  goes on top of the plane and the hydrophilic solvent goes below.

        We need to add some term to be consistent with the chosen Cs values.
        csi > csj means that csi is more hydrophilic than csj
   
    '''
    _instance = None
    def __new__(cls): 
        '''
        This function (SINGLETON pattern) has been created to force the class to be initialised only once and to be able to access the stored values, to read cs1,cs2,w1,w2 
        without modifying much the initial code and that the user can give them in the main.py as input.
        '''
        if cls._instance is None:
            cls._instance = super(interface_parameters, cls).__new__(cls)
            # Initial value must be positive!
            cls._instance._weight_interface = 0.1 #0.1
            cls._instance._last_calculated_weight = None
        return cls._instance 
        
    
    def calculate_weight_interface (self, cs_phase_1: float, cs_phase_2: float, ws_phase_1: float ,ws_phase_2: float):
        # Be consistent with th chosen Cs values
        if cs_phase_1 == cs_phase_2 and ws_phase_1 == ws_phase_2:
            mj_modified_parameters = 0. # No interface!
        elif  cs_phase_1 == cs_phase_2 and ws_phase_1 != ws_phase_2:
            mj_modified_parameters = -1
        elif cs_phase_1 > cs_phase_2:
            mj_modified_parameters = -1
        elif cs_phase_2 > cs_phase_1:
            mj_modified_parameters = 1
        else:
            mj_modified_parameters = 1
        
        self._last_calculated_weight = mj_modified_parameters*self._weight_interface
        print('Interface weight: ',self._last_calculated_weight)
        return  self._last_calculated_weight#mj_modified_parameters*self._weight_interface 
    


    # Used to access the value stored in other classes where the instance has not been initialised.
    @property
    def get_weight_interface(self):
        return  self._last_calculated_weight 
