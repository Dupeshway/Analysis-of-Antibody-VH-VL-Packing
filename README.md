README: VH/VL AB Packing Angle Analysis Project
Author: Yobi Livingstone

Essential: Execute each Section separately to generate dataset, dataframes and plots.
Reasoning: Some sections require copying contents into BASH or light formatting
in spreadsheet software such as Libre Calculation/ Excel/ Google Sheets. Within each section in the .py files, there are instructions on how to generate the data.
Sections:
Generating and cleaning the pdb packing angle dataset from AbDb (parsing_suite.py)
Generating an All dataset of packing angles
Extract Free vs Complex pdbs for Packing angles
Formatting All, Free and Complex for RStudio
Generating scripts to derive mean, range, sd for RStudio
Dataframes for Pandas and Seaborn
Flexibility grade: extracting PDBs by list in textfile
Extracting Residues
Generating scripts for data frames (df_pangles.py)
Plotting data frames: Mean, Median, Range, SD: (df_plots.py)
Melting data frames
Mean free vs comp across all sets
Mean of largest redundancies
Standard Deviations
Range: Free / Comp / All (density)
Range: Free / Comp / All (swarm/violin plot)
Range redundancy plot (all 4-6 and 7+)
Range density plot (all 3 in 1)
Range against redund count (all 3 in 1)
Standard deviations
SD of Density plot(all 3 in 1)
SD and Range Scatter plot
SD Range Free vs Complex
Archive
Flexibility gradients (df_flexibility.py)
Extracting residues
Grade 1 extraction
Grade 2 extraction
Grade 3 extraction
Grade R extraction
Plotting Hydrophobicity: (df_hydroplots.py)
Melting Dataframes: isolating columns and melting data
Total vs Flexibility: Plotting Hydrophobic values against flexibility gradients
L and H correlation: Plotting Light and Heavy chain hydrophobic values 
Weighted L and H correlation: Plotting light and heavy chains weighted by the probability of residues position being a contact residue 
Modules created for extracting, cleaning, analysing and plotting dataset
Extracting and cleaning antibody packing angles for data frames (cleaning_pdb)
Extracting residues from positions in .pdb file (residue_extract.py)
Opening, reading and writing functions (file_managment.py)
Configuring pathways and files (config.py)

 

