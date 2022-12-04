# Stack_detect
A python code for detecting stacking interaction in RNA structures
 The code works on Python3. To run the code in a terminal go through the following steps.

Download the whole code into a folder
Unzip the "utils" folder and keep it in main folder
Create a subfolder (eg PDBfiles) and store the target PDB structurs into it.
In the terminal run the code by >>>Python run_stacking.py PDBfile/ This will create files for each PDB structure in the PDBfiles folder.
To consolidate all the data into a single text files ">>>Python consolidate_stack.py PDBfile/ ". This will result in a  raw "output.txt"
To get the final output, run the code:  ">>>Python getMcAnnotate_notation.py output.txt"      This will porduce the final result as  "final_output.txt"
                                        
