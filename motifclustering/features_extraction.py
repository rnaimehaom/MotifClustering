import sys, os
import numpy as np
import pandas as pd 
import subprocess


from Bio.SeqUtils.ProtParam import ProteinAnalysis



class MotifProperties:
    """
    Extract features/measurments from motifs sequences. 
    Rely on Bio.SeqUtils.ProtParam package.
    see https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html. for 
    more informations.
    
    
    The computed features (columns of the output dataframe) are: 

    - 'id' (str) :                  the motif
    - 'seq_len' (int) :             the motif length 
    - 'mol_wt' (float) :            the total molecular weight
    - 'gravy' (float) :             (grand average of hydropathy) calculated by 
                                    adding the hydropathy value for each residue 
                                    and dividing by the length of the sequence 
                                    (Kyte and Doolittle; 1982) 
                                    (Bio.SeqUtils.ProtParam documentation).
                                    NB : probably redundant with 'mean_kd_hydro' 
                                    and one should be remooved in the futur.
    - 'isoelec_point' (float) :     calculate the isoelectric point. 
                                    (Bio.SeqUtils.ProtParam documentation)
    - 'charge_at_pH' (float) :      calculate the charge of a protein at given 
                                    pH (default : 7.5).
                                    (Bio.SeqUtils.ProtParam documentation) 
                                    "The intracellular pH of living cells is 
                                    strictly controlled in each compartment. 
                                    Under normal conditions, the cytoplasmic pH 
                                    (pHc) and the vacuolar pH (pHv) of typical 
                                    plant cells are maintained at slightly 
                                    alkaline (typically 7.5) and acidic 
                                    (typically 5.5) values, respectively 
                                    (https://doi.org/10.1007/978-3-7091-1254 \
                                    -0_4).
    - 'molar_ext_coef' (float) :    calculates the molar extinction coefficient 
                                    assuming cysteines (reduced) and cystines 
                                    residues (Cys-Cys-bond)
                                    (Bio.SeqUtils.ProtParam documentation)
    - 'tiny' (float) :              cumulative percentage of the smallest amino 
                                    acids (e.g 'A', 'C', 'G', 'S', 'T')
    - 'small' (float) :             cumulative percentage of the small amino 
                                    acids (e.g A', 'C', 'F', 'G', 'I', 'L', 'M', 
                                    'P', 'V', 'W', 'Y').
    - 'aliphatic' (float) :         cumulative percentage of the aliphatic amino 
                                    acids (e.g 'A', 'I', 'L', 'V')
    - 'aromatic' (float) :          cumulative percentage of the aromatic amino 
                                    acids (e.g 'F', 'H', 'W', 'Y')
    - 'non_polar' (float) :         cumulative percentage of the non-polar amino 
                                    acids (e.g 'A', 'C', 'F', 'G', 'I', 'L', 
                                    'M', 'P', 'V', 'W', 'Y')
    - 'polar' (float) :             cumulative percentage of the polar amino 
                                    acids (e.g 'D','E', 'H', 'K', 'N', 'Q', 'R', 
                                    'S', 'T', 'Z')
    - 'charged' (float) :           cumulative percentage of the charged amino 
                                    acids (e.g 'B', 'D', 'E', 'H', 'K', 'R', 
                                    'Z')
    - 'basic' (float) :             cumulative percentage of the basic amino 
                                    acids (e.g 'H', 'K', 'R')
    - 'acidic' (float) :            cumulative percentage of the acidic amino 
                                    acids (e.g 'B', 'D', 'E', 'Z')
    - 'helix' (float) :             cumulative percentage of amino acids which 
                                    tend to be in Helix (e.g. 'V', 'I', 'Y', 
                                    'F', 'W', 'L')
    - 'turn' (float) :              cumulative percentage of amino acids which 
                                    tend to be in turn (e.g 'N', 'P', 'G', 'S')
    - 'sheet' (float) :             cumulative percentage of amino acids which 
                                    tend to be in turn (e.g  'E', 'M', 'A', 'L')
    """
    
    def __init__(self):
        self.lst_prot_prop = []
        self.dict_prop = {
            'tiny':['A', 'C', 'G', 'S', 'T'],
            'small':['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
            'aliphatic':['A', 'I', 'L', 'V'],
            'aromatic':['F', 'H', 'W', 'Y'],
            'non_polar':['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
            'polar':['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z'],
            'charged':['B', 'D', 'E', 'H', 'K', 'R', 'Z'],
            'basic':['H', 'K', 'R'],
            'acidic':['B', 'D', 'E', 'Z'],
            'kyte_doolittle':{'A': 1.8,'C': 2.5,'D': -3.5,'E': -3.5,
                              'F': 2.8,'G': -0.4,'H': -3.2,'I': 4.5,
                              'K': -3.9,'L': 3.8,'M': 1.9,'N': -3.5,
                              'P': -1.6,'Q': -3.5,'R': -4.5,'S': -0.8,
                              'T': -0.7,'V': 4.2,'W': -0.9,'Y': -1.3},
            'rose':{'A': 0.74,'C': 0.91,'D': 0.62,'E': 0.62,'F': 0.88,
                    'G': 0.72,'H': 0.78,'I': 0.88,'K': 0.52,'L': 0.85,
                    'M': 0.85,'N': 0.63,'P': 0.64,'Q': 0.62,'R': 0.64,
                    'S': 0.66,'T': 0.70,'V': 0.86,'W': 0.85,'Y': 0.76}}
        
    
    def extract_properties(self, lst_motif, ph=7.5):
        """
        Iterate through a list of motif sequences to compute 
        several protein properties
        """
        
        self.ph = ph
        for motif in lst_motif:
            self.compute_properties(motif)
            
            
    def compute_properties(self, motif):
        """
        Compute properties from the sequence.
        The computed properties are the ones described in the class 
        description.
        """
        motif_analysis = ProteinAnalysis(motif)
        motif_len=len(motif)
        aa_percent = motif_analysis.get_amino_acids_percent()
        sec_struct_fraction = motif_analysis.secondary_structure_fraction()
        self.lst_prot_prop.append({
            'motif':motif,
            'seq_len':motif_len,
            'avg_mol_wt':motif_analysis.molecular_weight()/motif_len,
            'gravy':motif_analysis.gravy(),
            'isoelec_point':motif_analysis.isoelectric_point(),
            'charge_at_pH':motif_analysis.charge_at_pH(self.ph),
            'molar_ext_coef':motif_analysis.molar_extinction_coefficient()[1],
            'tiny':np.sum([aa_percent[aa] for aa in self.dict_prop['tiny'] if \
                aa in aa_percent]),
            'small':np.sum([aa_percent[aa] for aa in self.dict_prop['small'] \
                if aa in aa_percent]),
            'aliphatic':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['aliphatic'] if aa in aa_percent]),
            'aromatic':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['aromatic'] if aa in aa_percent]),
            'non_polar':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['non_polar'] if aa in aa_percent]),
            'polar':np.sum([aa_percent[aa] for aa in self.dict_prop['polar'] \
                if aa in aa_percent]),
            'charged':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['charged'] if aa in aa_percent]),
            'basic':np.sum([aa_percent[aa] for aa in self.dict_prop['basic'] \
                if aa in aa_percent]),
            'acidic':np.sum([aa_percent[aa] for aa in self.dict_prop['acidic'] \
                if aa in aa_percent]),
            'small':np.sum([aa_percent[aa] for aa in self.dict_prop['small'] \
                if aa in aa_percent]),
            'helix':sec_struct_fraction[0],
            'turn':sec_struct_fraction[1],
            'sheet':sec_struct_fraction[2]
        })

    def export_as_dataframe(self):
        """
        Return a pandas dataframe from a list of dict.
        """
        
        return pd.DataFrame(self.lst_prot_prop)