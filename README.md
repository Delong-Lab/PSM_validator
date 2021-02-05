# PSM_validator: instructions for use

### Introduction

PSM Validation with Internal Standards (P-VIS) is a workflow designed to systematically and objectively assess the validity of invidiual peptide-spectrum matches (PSM) generated by analysis of tandem mass spectrometry data. PSM_validator is a computer program that automates the data analysis portion of the P-VIS workflow. For detailed information on P-VIS and PSM_validator, see [Peptide-Spectrum Match Validation with Internal Standards (P-VIS): Internally-Controlled Validation of Mass Spectrometry-Based Peptide Identifications (Wiles TA et al, 2021, J Proteome Res)](https://doi.org/10.1021/acs.jproteome.0c00355).

PSM_validator code is written in Python 3.7.4. Both the PSM_validator source code and a standalone executable (generated using PyInstaller 3.5) are available in the [release](https://github.com/Delong-Lab/PSM_validator/releases) section of this repository. To execute the source code, you will need the msconvert tool (part of ProteoWizard) installed on your computer. You will also need a Python interpreter and the Matplotlib and SciPy libraries; these libraries are included in the Anaconda distribution of Python. If using the standalone PSM_validator executable, you will only need to download the executable and ProteoWizard.

### Using PSM_validator

- [Download](http://proteowizard.sourceforge.net/) ProteoWizard.
- [Download](https://github.com/Delong-Lab/PSM_validator/releases) the latest release of PSM_validator and unzip the files.
- For each PSM being validated, locate the raw data file for the biological sample and the raw data file for the validation sample.
- Place all raw data files in a single directory. The data files and PSM_validator do not have to be in the same directory. 
- `PSM_validator.exe` (or the source code) is associated with a directory titled `parameters`. `PSM_validator.exe` (or the source code files) and the `parameters` directory must be kept together within the same directory for the program to run. 
- The `parameters` directory contains five CSV files that allow for user input. The files are most easily manipulated in a spreadsheet editor such as Microsoft Excel. 
  - `amino_acids.csv`: A list of the 20 canonical amino acids (single-letter code) and their monoisotopic residue masses. The user can add custom amino acids to this list by providing a single-letter name and an exact mass, thus allowing PSM_validator to handle modified amino acids. 
  - `msconvert_settings.csv`: The user indicates whether data files are already in MS1 and MGF formats or need to be converted from raw data files using msconvert. The user MUST indicate the directory where `msconvert.exe` is located. This should be in the main `ProteoWizard` directory. The input file type for generating MS1 and MGF files should be entered; usually, the same file type will serve as the input for generating both MS1 and MGF files. The file type must be one compatible with msconvert, such as Thermo .raw or Agilent .d. The user can enter a list of filters just as they would be entered in the msconvert command-line interface. Using appropriate filters, such as removing low-intensity peaks from the data, can dramatically reduce the subsquent PSM_validator analysis time. Note that a space character must be kept at the beginning of the filters in msconvert_settings.csv. For details about msconvert filters, visit http://proteowizard.sourceforge.net/tools/filters.html.
  - `queries.csv`: Where the user enters the analyses to be performed in batch format. The user indicates the directory containing all raw data files. For each analysis, the user indicates the name of the biological file and the validation file and the putative peptide sequence (the query sequence) being validated. If the sequence contains an N- or C-terminal modification, this can be indicated by entering the mass shift. Multiple analyses are performed automatically in the order listed. Selecting the verbose option generates figures for analyses of internal standard peptides (ISPs).
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

PSM_validator is freely available under the [Creative Commons Attribution 4.0 International Public License](https://creativecommons.org/licenses/by/4.0/legalcode). To provide acknowledgment, please cite the following manuscript: [Peptide-Spectrum Match Validation with Internal Standards (P-VIS): Internally-Controlled Validation of Mass Spectrometry-Based Peptide Identifications (Wiles TA et al, 2021, J Proteome Res)](https://doi.org/10.1021/acs.jproteome.0c00355). 


