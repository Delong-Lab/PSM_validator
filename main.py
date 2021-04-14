version = "PSM_validator_v1p6"

##################################################################################################################################################################

from os import path, system, listdir, remove
from csv import reader, writer
from datetime import datetime
from auxiliary import time_format
from refinement import precursor_refine, spectrum_merge
from validation import validation

##################################################################################################################################################################

print()
print("###############################################################################################################################")
print("#     PSM_VALIDATOR     #     PSM_VALIDATOR     #     PSM_VALIDATOR     #     PSM_VALIDATOR     #     PSM_VALIDATOR     #      ")   
print("###############################################################################################################################")
print()

scriptpath = path.realpath(__file__)
scriptdir = path.dirname(scriptpath)

now = datetime.now()
now = time_format(now)
date, time = now[0], now[1]
timestamp = date + "_" + time  

settings = [] 
settings_file = open(scriptdir+"\\parameters\\settings.csv")
settings_contents=reader(settings_file, delimiter=",")
row_num = 0
for row in settings_contents:
    if 2 <= row_num <= 6 or 21 <= row_num:
        parameter = list(row)
        settings.append(parameter[2])
    elif row_num > 1:
        parameter = list(row)
        settings.append(float(parameter[2]))
    row_num = row_num + 1
settings_file.close()

convert = settings[0]
msconvert_directory = settings[1]
filetype = settings[2]
MS1_filters = settings[3]
MGF_filters = settings[4]
pre_mz_tol = settings[5]
pro_mz_tol = settings[6]  
RTtol = settings[13]
space_saver = settings[20]
verbose = settings[21]
   
queries_file = open(scriptdir+"\\parameters\\queries.csv")
queries_contents = reader(queries_file, delimiter = ",")
row_num = 0
for row in queries_contents:
    if row_num == 1:
        directory_row = list(row)
        directory = directory_row[1]
        
        if convert == "Y":
            
            files = listdir(directory)
            for i in range(0, len(files)):
                if ".mgf" in files[i] or ".ms1" in files[i]:
                    remove(directory + "\\" + files[i])
                    
            print("Converting/filtering data files...")
            print()
            MS1_settings = 'cd ' + msconvert_directory + '& msconvert ' + directory + '\*' + filetype + ' --ms1 -o ' + directory + ' ' + MS1_filters
            MGF_settings = 'cd ' + msconvert_directory + '& msconvert ' + directory + '\*' + filetype + ' --mgf -o ' + directory + ' ' + MGF_filters
            command1 = 'cmd /c "' + MS1_settings + '"'
            command2 = 'cmd /c "' + MGF_settings + '"'
            system(command1)
            system(command2)
            msconvert_settings = [[MS1_settings], [MGF_settings], ["Precursor refinement and MS2 spectrum merging performed"]] 
            
            print("Performing precursor refinement...")
            print()
            files = listdir(directory)
            for i in range(0, len(files)):
                if ".mgf" in files[i]:
                    sample_mgf = directory + "\\" + files[i]
                    sample_name = files[i][0:len(files[i])-4]
                    sample_ms1 = directory + "\\" + sample_name + ".ms1"
                    precursor_refine(sample_ms1, sample_mgf)
                    
            print("Merging MS2 spectra...")
            print()
            files = listdir(directory)
            for i in range(0, len(files)):
                if "_precursors_refined.mgf" in files[i]:
                    refined_mgf = directory + "\\" + files[i]
                    spectrum_merge(refined_mgf, RTtol, pre_mz_tol, pro_mz_tol)
                    
        else: 
            msconvert_settings = [["File conversion/filtering, precursor refinement, and MS2 spectrum merging not performed"]]

        if space_saver == "Y":
            files = listdir(directory)
            for i in range(0, len(files)):
                if "_precursors_refined.mgf" in files[i] or "_raw.mgf" in files[i]:
                    remove(directory + "\\" + files[i])
            
    if row_num > 2:   
        analysis = list(row)
        if analysis[0] != "":
            if row_num == 3:
                batch_results=open(directory + "\\PSM_validator_" + timestamp + "_results.csv","w",newline="")
                batch_results_writer=writer(batch_results) 
                batch_results_writer.writerow(["DIRECTORY:", directory])
                batch_results_writer.writerow(["bio sample", "syn sample", "N-term mass shift", "sequence", "C-term mass shift", "processing time (min:sec)", "PCC", "PCC percentile rank", "deviation from expected RT (minutes)", "RT percentile rank", "PCC outcome", "RT outcome", "WARNINGS"])
            biological, synthetic = analysis[0], analysis[1]
            sequence = analysis[3]
            N_term_shift, C_term_shift = float(analysis[2]), float(analysis[4])
            results = validation(scriptdir, sequence, N_term_shift, C_term_shift, directory, biological, synthetic, verbose, msconvert_settings, version)
            output = analysis + list(results)
            batch_results_writer.writerow(output)
    row_num = row_num + 1
queries_file.close()
batch_results.close()

print("###############################################################################################################################")
print("All of your analyses are complete. Press any key to exit.")
print("###############################################################################################################################")
input()
