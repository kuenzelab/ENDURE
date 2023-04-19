# Welcome to ENDURE! 💻🔬

**ENDURE is a user-friendly web application designed to help you analyze the energetic contributions of single and multi-variant protein designs. With ENDURE, you can quickly and easily evaluate the structural and thermodynamic changes that occur when mutations are introduced into your protein designs.**

## File Upload 📂

The first step in using ENDURE is to upload your protein models in PDB format. On the file upload page, simply select the PDB files you want to analyze and hit the "Upload" button. ENDURE will take care of the rest, running the necessary preprocessing steps to prepare your files for analysis.

## Interaction Analysis 🔍

This section provides an overview of the interactions between the residues in the uploaded protein structure. The user can select different types of interactions, such as salt bridges, sulfide bonds, and hydrogen bonds, and view their changes between a variant and wild-type structure.

## Residue Depth 📈

The residue depth page provides a visual representation of the change in residue depth that occurs when mutations are introduced. Residue depth is calculated using the Biopython library and is defined as the average distance (in angstroms) of the atoms in a residue from the solvent accessible surface. By analyzing the changes in residue depth, you can gain insights into how mutations affect the structural stability of your designs.

## Energy Heatmap 🔥

The energy heatmap page provides a visual representation of the changes in interaction energy that occur when mutations are introduced. The heatmap allows you to easily identify which residues are contributing the most to the changes in interaction energy, providing valuable insights into how to optimize your designs for stability and functionality.

## Example user workflow

In order to utilize the full capabilities of the ENDURE tool, it is important to understand the sequential nature of the tool and follow a specific workflow. A typical user would start by uploading their protein model in PDB format on the File Upload page. Once the model has been uploaded, the preprocessing steps including cleaning the PDB files, determining mutations, and calculating residue depth must be performed. Once these steps have been completed, the user can then proceed to the Interaction Analysis page, where the energetic contributions of single and multi-variant protein designs can be analyzed. The Residue Depth page provides a visual representation of the average distance of atoms in a residue from the solvent accessible surface, while the Energy Heatmap page displays the energy contribution of each residue in the protein model.

It is recommended for new users to follow this pipeline in order to effectively utilize the ENDURE tool:

1. Upload protein model in PDB format on the File Upload page.
2. Perform preprocessing steps including cleaning PDB files, determining mutations, and calculating residue depth.
3. Analyze the energetic contributions of single and multi-variant protein designs on the Interaction Analysis page.
4. View the Residue Depth and Energy Heatmap pages for additional insights into the protein model.

## Foreword 📖
We thank Paola Engelberger-Aliaga for the draft of the tool's logo! 🎨

With ENDURE, you have all the tools you need to make informed decisions about your protein designs. So why wait? Get started today and start exploring the exciting world of protein design! 🚀

We are constantly working to improve and expand the capabilities of ENDURE. We welcome suggestions for new functionalities and any feedback you may have about the tool. If you encounter any issues or would like to make a suggestion, please raise an issue in our GitHub repository. We are always looking for ways to make ENDURE the most helpful and efficient tool for analyzing protein designs.

felipeengelberger@gmail.com

Thank you for using ENDURE! 💻🔬
