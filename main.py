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


#in readme have to tell user to install proteowizard
#also use msconvert spectrum merging feature and see difference it makes; can doing this make results as good as for specmill mzxml?
#it's really important that precursor refine works right; check that it does

        
queries_file = open(scriptdir+"\\parameters\\queries.csv")
queries_contents = csv.reader(queries_file, delimiter = ",")
row_num = 0
for row in queries_contents:
    if row_num == 1:
        directory_row = list(row)
        directory = directory_row[1]
        print("Converting files...")
        print()
         
        msconvert_settings_file = open(scriptdir+"\\parameters\\msconvert_settings.csv")
        msconvert_settings_contents=csv.reader(msconvert_settings_file, delimiter=",")
        line = 0
        for row in msconvert_settings_contents:
            if line == 2:
                msconvert_directory_row = list(row)
                msconvert_directory = msconvert_directory_row[1]
            elif line == 5:
                MS1_filetype_row = list(row)
                MS1_filetype = MS1_filetype_row[1]
            elif line == 6:
                MS1_filters_row = list(row)
                MS1_filters = MS1_filters_row[1]
            elif line == 9:
                MGF_filetype_row = list(row)
                MGF_filetype = MGF_filetype_row[1]
            elif line == 10:
                MGF_filters_row = list(row)
                MGF_filters = MGF_filters_row[1]
            line = line + 1
        msconvert_settings_file.close()
        
        MS1_settings = 'cd ' + msconvert_directory + '& msconvert ' + directory + '\*' + MS1_filetype + ' --ms1 -o ' + directory + MS1_filters
        MGF_settings = 'cd ' + msconvert_directory + '& msconvert ' + directory + '\*' + MGF_filetype + ' --mgf -o ' + directory + MGF_filters
        command1 = "cmd /c " + MS1_settings
        command2 = "cmd /c " + MGF_settings
        os.system(command1)
        os.system(command2)
        msconvert_settings = [[MS1_settings], [MGF_settings]]      
    if row_num > 2:   
        analysis = list(row)
        if analysis[0] != "":
            if row_num == 3:
                batch_results=open(directory + "\\PSM_validator_" + timestamp + "_results.csv","w",newline="")
                batch_results_writer=csv.writer(batch_results) 
                batch_results_writer.writerow(["DIRECTORY:", directory])
                batch_results_writer.writerow(["bio sample", "syn sample", "N-term mass shift", "sequence", "C-term mass shift", "verbose", "processing time (min:sec)", "PCC", "PCC percentile rank", "deviation from expected RT (minutes)", "RT percentile rank", "PCC outcome", "RT outcome", "WARNINGS"])
            biological, synthetic = analysis[0], analysis[1]
            sequence = analysis[3]
            N_term_shift, C_term_shift = float(analysis[2]), float(analysis[4])
            verbose = analysis[5]
            results = validation(scriptdir, sequence, N_term_shift, C_term_shift, directory, biological, synthetic, verbose, msconvert_settings)
            output = analysis + list(results)
            batch_results_writer.writerow(output)
    row_num = row_num + 1
queries_file.close()
batch_results.close()

print("###############################################################################################################################")
print("All of your analyses are complete. Press any key to exit.")
print("###############################################################################################################################")
input()
