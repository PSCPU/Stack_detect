A python code for detecting stacking interaction in RNA structures The code works on Python3. To run the code in a terminal go through the following steps.

1.Download the whole code into a folder

2.Unzip the "utils" folder and keep it in main folder 

3.Create a subfolder (eg PDBfiles) and store the target PDB structurs into it. 

4.In the terminal run the code by >>>Python run_stacking.py PDBfile/ This will create files for each PDB structure in the PDBfiles folder. 

5.To consolidate all the data into a single text files ">>>Python consolidate_stack.py PDBfile/ ". This will result in a raw "output.txt" 

6.To get the final output, run the code: ">>>Python getMcAnnotate_notation.py output.txt" This will porduce the final result as "final_output.txt"
