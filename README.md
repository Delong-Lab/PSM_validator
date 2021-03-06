# PSM_validator: instructions for use

### Introduction

PSM Validation with Internal Standards (P-VIS) is a workflow designed to systematically and objectively assess the validity of invidiual peptide-spectrum matches (PSM) generated by analysis of tandem mass spectrometry data. PSM_validator is a computer program that automates the data analysis portion of the P-VIS workflow. For detailed information on P-VIS and PSM_validator, see [cite paper].

PSM_validator code is written in Python 3.7.4. Both the PSM_validator source code and a standalone executable (generated using PyInstaller 3.5) are available in the [release](https://github.com/Delong-Lab/PSM_validator/releases) section of this repository. Dependencies for executing the source code are the Matplotlib and SciPy libraries, which are included in the Anaconda distribution of Python.

### Using PSM_validator

- [Download](https://github.com/Delong-Lab/PSM_validator/releases) the latest release of PSM_validator and unzip the files.
- For each PSM being validated, locate the raw data file for the biological sample and the raw data file for the validation sample.
- Using MSconvert (part of the [ProteoWizard](http://proteowizard.sourceforge.net/download.html) toolkit), convert each raw mass spectrometry data file to both a Mascot Generic Format (MGF) file and an MS1 file.
  - The MGF file and the MS1 file for a given sample must have the same name. 
  - All four files associated with a single analysis--the MGF file and the MS1 file for the biological sample and the MGF file and the MS1 file for the validation sample--must be in the same directory.
- `PSM_validator.exe` (or the source code) is associated with a directory titled `parameters`. `PSM_validator.exe` (or the source code files) and the `parameters` directory must be kept together within the same directory for the program to run. 
- The `parameters` directory contains four CSV files that allow for user input. The files are most easily manipulated in a spreadsheet editor such as Microsoft Excel. 
  - `amino_acids.csv`: A list of the 20 canonical amino acids (single-letter code) and their monoisotopic residue masses. The user can add custom amino acids to this list by providing a single-letter name and an exact mass, thus allowing PSM_validator to handle modified amino acids. 
  - `queries.csv`: Where the user enters the analyses to be performed in batch format. For each analysis, the user indicates the name of the biological file and the validation file, the directory containing the files, and the putative peptide sequence (the query sequence) being validated. If the sequence contains an N- or C-terminal modification, this can be indicated by entering the mass shift. Multiple analyses are performed automatically in the order listed. 
  - `settings.csv`: Provides access to all user-adjustable algorithm parameters, such as mass tolerances, along with a description of each and information on acceptable values. PSM_validator is capable of processing data generated by collision-induced dissociation (CID), which is dominated by b- and y-ions, or electron transfer dissociation (ETD), which is dominated by c- and z-ions. 
  - `standards.csv`: Contains the amino acid sequence for each of the internal standard peptides (ISPs). By default, the [PROCAL](https://shop.jpt.com/59-Proteomics-Peptides-Kits/140-Standardization-Kits/) sequences are listed. 
- After the files in `parameters` have been updated, saved, and closed, the user clicks on the `PSM_validator.exe` icon (or runs `main.py` using a Python interpreter) to begin the analyses.

### PSM_validator output

- Execution of PSM_validator opens a console window in which real-time results are displayed, allowing the user to monitor progress and results. 
- Output files for all analyses listed in `queries.csv` are saved in the directory containing the MGF and MS1 data files used for the first analysis. 
  - For each individual analysis, a time-stamped folder named after the query sequence is generated. This folder contains a detailed results summary file, which contains most of the information displayed in the real-time output, and two folders: `Figures` and `Tables`. 
    - `Figures`: Contains six automatically-generated, publication-ready plots. These plots provide an intuitive representation of the results and allow the user to manually check the veracity of the assessment and evaluate if appropriate settings were used. 
    - `Tables`: Contains detailed scoring results and other details from the analysis in tabular form.
  - After all analyses are complete, a time-stamped, high-level batch results summary file is generated containing the final results for each analysis. 

### Test data and results

Benchmarking data and associated PSM_validator result files have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD018370.

### License

PSM_validator is freely available under the [Creative Commons Attribution 4.0 International Public License](https://creativecommons.org/licenses/by/4.0/legalcode). To provide acknowledgment, please cite the following manuscript: https://doi.org/10.1021/acs.jproteome.0c00355. 


