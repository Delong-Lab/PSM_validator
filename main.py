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
    if row_num == 16:
        parameter = list(row)
        settings.append(parameter[2])
    elif row_num > 1:
        parameter = list(row)
        settings.append(float(parameter[2]))
    row_num = row_num + 1
settings_file.close()

pre_mz_tol = settings[0]
pro_mz_tol = settings[1]
abund_thresh = settings[2]
PCC_abund_thresh = settings[3]
min_score = settings[4]
min_weighted_score = settings[5]
min_pairs_PCC = settings[6]    
min_PCC = settings[7]    
RTtol = settings[8]
min_RT = settings[9]
max_RT = settings[10]
manual_RTdev_thresh = settings[11]
min_intstd = settings[12]
percentile_thresh = settings[13]
ion_type = settings[14]
       
queries_file = open(scriptdir+"\\parameters\\queries.csv")
queries_contents = reader(queries_file, delimiter = ",")
row_num = 0
for row in queries_contents:
    if row_num == 1:
        directory_row = list(row)
        directory = directory_row[1]
         
        msconvert_settings_file = open(scriptdir+"\\parameters\\msconvert_settings.csv")
        msconvert_settings_contents=reader(msconvert_settings_file, delimiter=",")
        line = 0
        for row in msconvert_settings_contents:
            if line == 2:
                msconvert_row = list(row)
                convert = msconvert_row[1]                   
            elif line == 4:
                msconvert_directory_row = list(row)
                msconvert_directory = msconvert_directory_row[1]
            elif line == 7:
                MS1_filetype_row = list(row)
                MS1_filetype = MS1_filetype_row[1]
            elif line == 8:
                MS1_filters_row = list(row)
                MS1_filters = MS1_filters_row[1]
            elif line == 11:
                MGF_filetype_row = list(row)
                MGF_filetype = MGF_filetype_row[1]
            elif line == 12:
                MGF_filters_row = list(row)
                MGF_filters = MGF_filters_row[1]
            line = line + 1
        msconvert_settings_file.close()
        
        if convert == "Y":
            
            files = listdir(directory)
            for i in range(0, len(files)):
                if ".mgf" in files[i] or ".ms1" in files[i]:
                    remove(directory + "\\" + files[i])
                    
            print("Converting/filtering data files...")
            print()
            MS1_settings = 'cd ' + msconvert_directory + '& msconvert ' + directory + '\*' + MS1_filetype + ' --ms1 -o ' + directory + MS1_filters
            MGF_settings = 'cd ' + msconvert_directory + '& msconvert ' + directory + '\*' + MGF_filetype + ' --mgf -o ' + directory + MGF_filters
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
            
    if row_num > 2:   
        analysis = list(row)
        if analysis[0] != "":
            if row_num == 3:
                batch_results=open(directory + "\\PSM_validator_" + timestamp + "_results.csv","w",newline="")
                batch_results_writer=writer(batch_results) 
                batch_results_writer.writerow(["DIRECTORY:", directory])
                batch_results_writer.writerow(["bio sample", "syn sample", "N-term mass shift", "sequence", "C-term mass shift", "verbose", "processing time (min:sec)", "PCC", "PCC percentile rank", "deviation from expected RT (minutes)", "RT percentile rank", "PCC outcome", "RT outcome", "WARNINGS"])
            biological, synthetic = analysis[0], analysis[1]
            sequence = analysis[3]
            N_term_shift, C_term_shift = float(analysis[2]), float(analysis[4])
            verbose = analysis[5]
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
