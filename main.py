import os, csv
from datetime import datetime
from auxiliary import time_format
from validation import validation

scriptpath = os.path.realpath(__file__)
scriptdir = os.path.dirname(scriptpath)

##################################################################################################################################################################

print()
print("###############################################################################################################################")
print("#     PSM_VALIDATOR     #     PSM_VALIDATOR     #     PSM_VALIDATOR     #     PSM_VALIDATOR     #     PSM_VALIDATOR     #      ")   
print("###############################################################################################################################")
print()

now = datetime.now()
now = time_format(now)
date, time = now[0], now[1]
timestamp = date + "_" + time  
        
queries_file = open(scriptdir+"\\parameters\\queries.csv")
queries_contents = csv.reader(queries_file, delimiter = ",")
row_num = 0
for row in queries_contents:
    if row_num > 1:   
        analysis = list(row)
        directory = analysis[0]
        if row_num == 2:
            batch_results=open(directory + "\\PSM_validator_" + timestamp + "_results.csv","w",newline="")
            batch_results_writer=csv.writer(batch_results) 
            batch_results_writer.writerow(["directory", "bio sample", "syn sample", "N-term mass shift", "sequence", "C-term mass shift", "verbose", "processing time (min:sec)", "PCC", "PCC percentile rank", "deviation from expected RT (minutes)", "RT percentile rank", "PCC outcome", "RT outcome", "WARNINGS"])
        biological, synthetic = analysis[1], analysis[2]
        sequence = analysis[4]
        N_term_shift, C_term_shift = float(analysis[3]), float(analysis[5])
        verbose = analysis[6]
        results = validation(scriptdir, sequence, N_term_shift, C_term_shift, directory, biological, synthetic, verbose)
        output = analysis + list(results)
        batch_results_writer.writerow(output)
    row_num = row_num + 1
queries_file.close()
batch_results.close()

print("###############################################################################################################################")
print("All of your analyses are complete. Press any key to exit.")
print("###############################################################################################################################")
input()
