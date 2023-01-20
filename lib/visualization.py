from typing import List, Dict, Iterable
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
import py3Dmol
from stmol import showmol
from io import StringIO
import streamlit as st


class WebStructure:
    """
    A PDB structure to be used in a py3Dmol Web Viewer
    """
    def __init__(self, file_name: str):
        """
        Initialize the PDB Structure
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        """
        self.parser = PDBParser()
        self.parser.QUIET = True
        self.pdb_file: StringIO =\
            st.session_state['File Upload'][f'pdb_{file_name}_clean']
        self.structure = self.__load_structure()
        self.text = self.__load_text()

    def __load_structure(self) -> Structure:
        """
        Parse PDB file using Biopython
        :return:
            Biopython structure object
        """
        structure = self.parser.get_structure(0, self.pdb_file)
        self.pdb_file.seek(0)
        return structure

    def __load_text(self) -> str:
        """
        Read Raw text of PDB File
        :return:
            A String of the PDB Contents
        """
        text = self.pdb_file.read()
        self.pdb_file.seek(0)
        return text


class WebViewer:
    """
    A py3Dmol viewer object to display PDB files in the browser
    """
    def __init__(self, width: int = 800, height: int = 500):
        """
        Create the viewer object
        :param width:
            The width of the viewer in pixels
        :param height:
            The height of the viewer in pixels
        """
        self.view = py3Dmol.view(width=width, height=height)
        self.structs: Dict[str, WebStructure] = {}
        self.id_struct: Dict[str, int] = {}
        self.width = width
        self.height = height

    def add_model(self, file_name: str) -> None:
        """
        Add a structure to the viewer
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :return:
            None
        """
        structure = WebStructure(file_name)
        self.structs[file_name] = structure
        self.id_struct[file_name] =\
            0 if not len(self.id_struct) else max(self.id_struct.values()) + 1
        self.view.addModel(structure.text, file_name)
        self.view.zoomTo()

    def show_cartoon(self, file_name: str, color: str) -> None:
        """
        Show the entire structure as a cartoon
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param color:
            The color of the cartoon
        :return:
            None
        """
        self.__set_style(
            {
                'model': self.id_struct[file_name]
            },
            {
                'cartoon': {
                    'color': color
                }
            }
        )

    def color_cartoon(self, file_name: str, resi: int, color: str) -> None:
        """
        Set a particular location of the cartoon to a specific color
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param resi:
            The residue location to be colored
        :param color:
            The desired color
        :return:
            None
        """
        self.__set_style(
            {
                'model': self.id_struct[file_name],
                'resi': resi
            },
            {
                'cartoon': {
                    'color': color
                }
            }
        )

    def __set_style(self, criteria: dict, style: dict) -> None:
        """
        Wrapper to the javascript method
        :param criteria:
            Selection criteria to control how the style is applied
        :param style:
            The visual style that is being applied
        :return:
            None
        """
        self.view.setStyle(criteria, style)

    def show_sc(
        self,
        file_name: str,
        resi: Iterable[int],
        color_stick: str,
        color_cartoon: str
    ) -> None:
        """
        Show the side chains of specific residues
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param resi:
            A list of target residue positions
        :param color_stick:
            The color of the side chains
        :param color_cartoon:
            The color of the cartoon
        :return:
            None
        """
        self.__set_style(
            {
                'model': self.id_struct[file_name],
                'resi': list(resi),
                'not': {
                    'atom': ['N', 'CA', 'C', 'O']
                }
            },
            {
                'stick': {
                    'color': color_stick
                }
            }
        )
        self.__set_style(
            {
                'model': self.id_struct[file_name],
                'resi': list(resi),
                'or': [
                    {
                        'atom': 'CA'
                    },
                    {
                        'atom': 'N',
                        'resn': 'PRO'
                    }
                ],
            },
            {
                'stick': {
                    'color': color_stick
                },
                'cartoon': {
                    'color': color_cartoon
                }
            }
        )

    def set_background(self, color: str) -> None:
        """
        Sets the Background color of the viewer
        :param color:
            The background color
        :return:
            None
        """
        self.view.setBackgroundColor(color)

    def center(self, file_name: str, resi: List[int]) -> None:
        """
        Center the viewer at a particular location
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param resi:
            The selection to be centered around
        :return:
            None
        """
        self.view.center(
            {
                'model': self.id_struct[file_name],
                'resi': resi
            }
        )

    def show_hydrogen(
        self,
        file_name: str,
        resi: Iterable[int],
        color_stick: str
    ) -> None:
        """
        An unfinished attempt to show the hydrogen atoms of a side chain
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param resi:
            The residue locations to show hydrogen atoms
        :param color_stick:
            The color of the hydrogen atom sticks
        :return:
            None
        """
        self.__set_style(
            {
                'model': self.id_struct[file_name],
                'resi': list(resi),
                'elem': 'H'
            },
            {
                'stick': {
                    'color': color_stick
                }
            }
        )

    def label_resi(self, file_name: str, resi: List[int]) -> None:
        """
        Show the labels of a residue
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param resi:
            The residues to be labeled
        :return:
            None
        """
        self.view.addResLabels(
            {
                'model': self.id_struct[file_name],
                'resi': resi
            },
            {
                'font': 'Arial',
                'fontColor': 'black',
                'showBackground': False
            }
        )

    def show(self) -> None:
        """
        Use the Streamlit showmol function to embed the py3Dmol viewer object
        in the streamlit body
        :return:
            None
        """
        showmol(self.view, width=self.width, height=self.height)
