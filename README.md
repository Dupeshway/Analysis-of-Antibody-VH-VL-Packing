README: VH/VL AB Packing Angle Analysis Project
Author: Yobi Livingstone

Essential: Execute each Section separately to generate dataset, dataframes and plots.
Reasoning: Some sections require copying contents into BASH or light formatting
in spreadsheet software such as Libre Calculation/ Excel/ Google Sheets. Within each section in the .py files, there are instructions on how to generate the data.


Sections:

1 Generating and cleaning the pdb packing angle dataset from AbDb (parsing_suite.py)
1.1 	Generating an All dataset of packing angles
1.2 	Extract Free vs Complex pdbs for Packing angles
1.3 	Formatting All, Free and Complex for RStudio
1.4 	Generating scripts to derive mean, range, sd for RStudio
1.5 	Dataframes for Pandas and Seaborn
1.6 	Flexibility grade: extracting PDBs by list in textfile
1.7 	Extracting Residues
1.8 	Generating scripts for data frames (df_pangles.py)

2. Plotting data frames: Mean, Median, Range, SD: (df_plots.py)
2.1 	Melting data frames
2.2	Mean free vs comp across all sets
2.3	Mean of largest redundancies
2.4	Standard Deviations
2.5	Range: Free / Comp / All (density)
2.6	Range: Free / Comp / All (swarm/violin plot)
2.7	Range redundancy plot (all 4-6 and 7+)
2.8	Range density plot (all 3 in 1)
2.9	Range against redund count (all 3 in 1)
2.10	Standard deviations
2.11	SD of Density plot(all 3 in 1)
2.12	SD and Range Scatter plot
2.13	SD Range Free vs Complex
2.14	Archive

3. Flexibility gradients (df_flexibility.py)
3.1	Extracting residues
3.2	Grade 1 extraction
3.3	Grade 2 extraction
3.4	Grade 3 extraction
3.5	Grade R extraction

4. Plotting Hydrophobicity: (df_hydroplots.py)
4.1	Melting Dataframes: isolating columns and melting data
4.2	Total vs Flexibility: Plotting Hydrophobic values against flexibility gradients
4.3	L and H correlation: Plotting Light and Heavy chain hydrophobic values 
4.4	Weighted L and H correlation: Plotting light and heavy chains weighted by the probability of residues position being a contact residue 

5. Modules created for extracting, cleaning, analysing and plotting dataset
5.1	Extracting and cleaning antibody packing angles for data frames (cleaning_pdb)
5.2	Extracting residues from positions in .pdb file (residue_extract.py)
5.3	Opening, reading and writing functions (file_managment.py)
5.4	Configuring pathways and files (config.py)

