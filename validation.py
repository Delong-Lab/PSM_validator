from os import makedirs
from csv import reader, writer
from copy import deepcopy
from math import atanh
from datetime import datetime
from auxiliary import time_format, timer
from scoring_prep import sequence_mass_calculator, LR_ion_calculator, pre_mz_filter
from scoring import LR_score
from retention_time import MS1_RTfilter, extraction, RT_prediction
from PCC_calculator import PCC_calculator
from plotting import swarm_plot, regression_plot, mirror_plot, EIC
from statistical_analysis import PCC_percentile_rank, RT_percentile_rank

##################################################################################################################################################################

def validation(scriptdir, sequence, N_term_shift, C_term_shift, directory, biological, synthetic, verbose, msconvert_settings, version):

    now = datetime.now()
    now = time_format(now)
    date, time = now[0], now[1]
    timestamp = date + "_" + time

    hyphen=sequence.find("-")
    
    sequence_formatted = sequence
    if N_term_shift != 0:
        sequence_formatted = "("+str(round(N_term_shift,2))+")"+sequence_formatted 
    if C_term_shift != 0:
        sequence_formatted = sequence_formatted+"("+str(round(C_term_shift,2))+")"
    
    bio_MGFpath = directory + "\\" + biological + ".mgf"    
    bio_MS1path = directory + "\\" + biological + ".ms1"  
    syn_MGFpath = directory + "\\" + synthetic + ".mgf"
    syn_MS1path = directory + "\\" + synthetic + ".ms1" 
    
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
    
    pre_mz_tol = settings[5]
    pro_mz_tol = settings[6]
    abund_thresh = settings[7]
    PCC_abund_thresh = settings[8]
    min_score = settings[9]
    min_weighted_score = settings[10]
    min_pairs_PCC = settings[11]    
    min_PCC = settings[12]    
    RTtol = settings[13]
    min_RT = settings[14]
    max_RT = settings[15]
    manual_RTdev_thresh = settings[16]
    min_intstd = settings[17]
    percentile_thresh = settings[18]
    ion_type = settings[19]
    
    amino_acids = {}
    amino_acids_file = open(scriptdir+"\\parameters\\amino_acids.csv")
    amino_acids_contents = reader(amino_acids_file, delimiter = ",")
    row_num = 0
    for row in amino_acids_contents:
        if row_num > 1:
            amino_acid = list(row)
            amino_acids[amino_acid[0]] = float(amino_acid[1])
        row_num = row_num + 1
    amino_acids_file.close()
    
    targets, N_term_shift_list, C_term_shift_list = [], [], []
    standards_file = open(scriptdir+"\\parameters\\standards.csv")
    standards_contents=reader(standards_file, delimiter=",")
    row_num = 0
    for row in standards_contents:
        if row_num > 1:
            standard = list(row)
            standard_formatted = [int(standards_contents.line_num), standard[1], sequence_mass_calculator(standard[1], amino_acids, float(standard[0]), float(standard[2])), "charge"] #line number and "charge" are placeholders to my knowledge
            targets.append(standard_formatted)
            N_term_shift_list.append(float(standard[0]))
            C_term_shift_list.append(float(standard[2]))
        row_num = row_num + 1
    standards_file.close()
    
    targets.append(["putative",sequence])
    
    #Calculate neutral precursor mass
    
    sequence_mass = sequence_mass_calculator(sequence, amino_acids, N_term_shift, C_term_shift)
    targets[(len(targets)-1)].append(sequence_mass)
    targets[(len(targets)-1)].append("charge") #get rid of this once you update dataframes or leave it then have the search write it
    
    makedirs(directory+"\\"+timestamp+"_"+sequence_formatted)
    out_dir=directory+"\\"+timestamp+"_"+sequence_formatted  
    
    results=open(out_dir + "\\" + sequence_formatted + "_" + timestamp + "_results.csv","w",newline="")
    results_writer=writer(results) 
    results_writer.writerow(["date (yymmdd):", date, "time (hhmmss):", time])   
    results_writer.writerow([])
    results_writer.writerow(["MSCONVERT SETTINGS"])
    results_writer.writerow([])
    results_writer.writerows(msconvert_settings)
    results_writer.writerow([])
    results_writer.writerow(["PSM_VALIDATOR SETTINGS"])
    results_writer.writerow([])
    results_writer.writerow(["N-terminus mass shift:", N_term_shift, "C-terminus mass shift:", C_term_shift])
    results_writer.writerow([])
    results_writer.writerow(["sequence:", sequence_formatted])
    results_writer.writerow([])
    results_writer.writerow(["neutral mass:", sequence_mass])
    results_writer.writerow(["", "+1", "+2", "+3", "+4"])
    results_writer.writerow(["m/z", (sequence_mass + 1*1.007276435)/1, (sequence_mass + 2*1.007276435)/2, (sequence_mass + 3*1.007276435)/3, (sequence_mass + 4*1.007276435)/4])
    results_writer.writerow(["Precursor mass tolerance (+/- ppm):", pre_mz_tol])
    results_writer.writerow(["Product ion mass tolerance (+/- ppm):", pro_mz_tol])
    results_writer.writerow(["Product ion abundance threshold (noise threshold for scoring and PCC calculation):", abund_thresh])
    results_writer.writerow(["Product ion abundance threshold for PCC calculation (% of max):", PCC_abund_thresh])
    results_writer.writerow(["Minimum allowable % backbone coverage:", min_score])
    results_writer.writerow(["Minimum allowable % "+ion_type+" intensity:", min_weighted_score])
    results_writer.writerow(["Minimum number of pairs for PCC calculation:", min_pairs_PCC])
    results_writer.writerow(["Minimum allowable PCC:", min_PCC])
    results_writer.writerow(["Window size for precursor extraction during RT determination (+/- minutes):", RTtol])
    results_writer.writerow(["Lower bound for RT analysis (minutes; based on biological sample run)", min_RT])
    results_writer.writerow(["Upper bound for RT analysis (minutes; based on biological sample run)", max_RT])
    results_writer.writerow(["Manual threshold for deviation from predicted RT (+/- minutes)", manual_RTdev_thresh])
    results_writer.writerow(["Minimum number of internal standards for percentile rank calculation:", min_intstd])
    results_writer.writerow(["Minimum allowable percentile (%):", percentile_thresh])
    results_writer.writerow(["Fragment ion type ('"+ion_type+"' for CID or 'c/z' for ETD):", ion_type])
    results_writer.writerow([])
    results_writer.writerow(["directory:", directory, "biological sample:", biological, "synthetic sample:", synthetic])
    results_writer.writerow([])
    
    #Start timer.
    
    start_time=datetime.now()
    print("###############################################################################################################################")
    print("ANALYSIS STARTED FOR", sequence_formatted+"...")
    print("###############################################################################################################################")
    print()
        
    #Calculate singly- and doubly-charged b- and y-ion masses and assign to L or R peptide if sequence is hybrid/spliced
    
    for i in range(0,len(targets)):
        if i == len(targets) - 1:
            ions=LR_ion_calculator(targets[i][1],amino_acids, N_term_shift, C_term_shift, ion_type)
        else:
            ions=LR_ion_calculator(targets[i][1],amino_acids, N_term_shift_list[i], C_term_shift_list[i], ion_type)        
        targets[i].append(ions[0])
        targets[i].append(ions[1])
        targets[i].append(ions[2])
        targets[i].append(ions[3])
    
        #Add a score of 0 as a final entry in each row of the L_ions and R_ions lists. Also add four zeros for obs m/z and four zeros for obs abundance.
    
        for j in range (len(targets[i][5])):
            targets[i][5][j]=targets[i][5][j]+[0,"","","","","","","","","",""]
        if targets[i][1].find("-") != -1:
            for k in range (len(targets[i][7])):
                targets[i][7][k]=targets[i][7][k]+[0,"","","","","","","","","",""]
    
    #Find and score spectra with the correct precursor m/z in the bio file.
            
    bio_MGF=open(bio_MGFpath)
    bio_contents=reader(bio_MGF, delimiter=" ")
    
    bio_hit_list=pre_mz_filter(targets, pre_mz_tol, bio_contents)
    
    bio_MGF.close()
        
    syn_MGF=open(syn_MGFpath)
    syn_contents=reader(syn_MGF, delimiter=" ")
    
    syn_hit_list=pre_mz_filter(targets, pre_mz_tol, syn_contents)
    
    syn_MGF.close()
    
    if len(bio_hit_list[len(bio_hit_list)-1])==8 and len(syn_hit_list[len(syn_hit_list)-1])==8:
        print("     ********************** Both the bio and syn run lack any spectra with correct precursor mass. ***********************")
        print()
        results_writer.writerow(["BOTH THE BIO AND SYN RUN LACK ANY SPECTRA WITH CORRECT PRECURSOR MASS"])
        results.close()
        processing_time = timer(datetime.now(), start_time)
        return(processing_time, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "BOTH THE BIO AND SYN RUN LACK ANY SPECTRA WITH CORRECT PRECURSOR MASS")    
    elif len(bio_hit_list[len(bio_hit_list)-1])==8:
        print("     ******************************* Bio run lacks any spectra with correct precursor mass. ******************************")
        print()
        results_writer.writerow(["BIO RUN LACKS ANY SPECTRA WITH CORRECT PRECURSOR MASS"])
        results.close()
        processing_time = timer(datetime.now(), start_time)
        return(processing_time, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "BIO RUN LACKS ANY SPECTRA WITH CORRECT PRECURSOR MASS")
    elif len(syn_hit_list[len(syn_hit_list)-1])==8:
        print("     ******************************* Syn run lacks any spectra with correct precursor mass. ******************************")
        print()
        results_writer.writerow(["SYN RUN LACKS ANY SPECTRA WITH CORRECT PRECURSOR MASS"])
        results.close()
        processing_time = timer(datetime.now(), start_time)
        return(processing_time, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "SYN RUN LACKS ANY SPECTRA WITH CORRECT PRECURSOR MASS")
    
    makedirs(directory+"\\"+timestamp+"_"+sequence_formatted+"\\Tables") 
    makedirs(directory+"\\"+timestamp+"_"+sequence_formatted+"\\Figures")
    makedirs(directory+"\\"+timestamp+"_"+sequence_formatted+"\\Figures\\normality") 
    if verbose == "Y":
        makedirs(directory+"\\"+timestamp+"_"+sequence_formatted+"\\Figures\\ISPs")         
    
    top_bio_hit=open(str(out_dir + "\\Tables\\" + sequence_formatted + "_bio_best_hit.mgf"),"w",newline="")
    top_bio_hit_writer=writer(top_bio_hit,delimiter=" ")
    bio_hits=open(str(out_dir + "\\Tables\\" + sequence_formatted + "_bio_hits.csv"),"w",newline="")
    bio_hits_writer=writer(bio_hits)
    
    print("SCORING SPECTRA FOR PEPTIDE OF INTEREST IN BIOLOGICAL SAMPLE...")
    results_writer.writerow(["QUERY RESULTS FOR BIOLOGICAL SAMPLE..."])
        
    top_bio_scans, bio_charge_list, bio_L_score, bio_R_score, bio_query_details = LR_score(abund_thresh,pro_mz_tol,min_score,min_weighted_score,bio_hit_list,bio_hits_writer,top_bio_hit_writer, "N", ion_type)
    
    top_bio_hit.close()
    bio_hits.close()
    
    query_bio_score = top_bio_scans[len(top_bio_scans)-1][1]
    query_bio_weighted_score= top_bio_scans[len(top_bio_scans)-1][2]
    
    results_writer.writerow([])
    if hyphen != -1:
        results_writer.writerow(["L pep backbone coverage", "R pep backbone coverage", "% backbone coverage", "% "+ion_type+" intensity", "charge"])
    else:
        results_writer.writerow(["% backbone coverage", "% "+ion_type+" intensity", "charge"])    
    for i in range(0, len(bio_query_details)):
        results_writer.writerow(bio_query_details[i])
    results_writer.writerow([])            
    
    if hyphen != -1:
        results_writer.writerow(["winning % backbone coverage:", round(query_bio_score, 1), "L pep backbone coverage:", "'"+bio_L_score, "R pep backbone coverage:", "'"+bio_R_score, "% "+ion_type+" intensity:", round(query_bio_weighted_score, 1), "charge:", bio_charge_list[len(bio_charge_list)-1]])
    else:
        results_writer.writerow(["winning % backbone coverage:", round(query_bio_score, 1), "% "+ion_type+" intensity:", round(query_bio_weighted_score, 1), "charge:", bio_charge_list[len(bio_charge_list)-1]])    
    results_writer.writerow([]) 
        
    #Find and score spectra with the correct precursor m/z in the syn file.
    
    top_syn_hit=open(str(out_dir + "\\Tables\\" + sequence_formatted + "_syn_best_hit.mgf"),"w",newline="")
    top_syn_hit_writer=writer(top_syn_hit,delimiter=" ")
    syn_hits=open(str(out_dir + "\\Tables\\" + sequence_formatted + "_syn_hits.csv"),"w",newline="")
    syn_hits_writer=writer(syn_hits)
    print("###############################################################################################################################")
    print()
    print("SCORING SPECTRA FOR PEPTIDE OF INTEREST IN SYNTHETIC SAMPLE...")
    results_writer.writerow(["QUERY RESULTS FOR SYNTHETIC SAMPLE..."])
    results_writer.writerow([])
    results_writer.writerow(["NOTE:", "Algorithm chooses winning syn spectrum only from those spectra with a corresponding precursor charge that matches"])
    results_writer.writerow(["", "the precursor charge for the winning spectrum chosen from the biological sample"])
    results_writer.writerow([])
    
    top_syn_scans, syn_charge_list, syn_L_score, syn_R_score, syn_query_details = LR_score(abund_thresh,pro_mz_tol,min_score,min_weighted_score,syn_hit_list,syn_hits_writer,top_syn_hit_writer, bio_charge_list, ion_type)
    
    top_syn_hit.close()
    syn_hits.close()
    
    query_syn_score= top_syn_scans[len(top_syn_scans)-1][1]
    query_syn_weighted_score= top_syn_scans[len(top_syn_scans)-1][2]

    results_writer.writerow([])
    if hyphen != -1:
        results_writer.writerow(["L pep backbone coverage", "R pep backbone coverage", "% backbone coverage", "% "+ion_type+" intensity", "charge"])
    else:
        results_writer.writerow(["% backbone coverage", "% "+ion_type+" intensity", "charge"])    
    for i in range(0, len(syn_query_details)):
        results_writer.writerow(syn_query_details[i])
    results_writer.writerow([])
    
    if hyphen != -1:
        results_writer.writerow(["winning % backbone coverage:", round(query_syn_score, 1), "L pep backbone coverage:", "'"+syn_L_score, "R pep backbone coverage:", "'"+syn_R_score, "% "+ion_type+" intensity:", round(query_syn_weighted_score, 1), "charge:", syn_charge_list[len(syn_charge_list)-1]])
    else:
        results_writer.writerow(["winning % backbone coverage:", round(query_syn_score, 1), "% "+ion_type+" intensity:", round(query_syn_weighted_score, 1), "charge:", syn_charge_list[len(syn_charge_list)-1]])    
    results_writer.writerow([])

    print("###############################################################################################################################",end="\n\n")
    print("CALCULATING PEARSON CORRELATION COEFFICIENTS (PCC) FOR INTERNAL STANDARDS...", end="\n\n")    
    print("[  sequence  ]            [bio % backbone coverage]  [bio % "+ion_type+" intensity]  [Pearson r]  [Fisher's z]  [# of pairs]  [z]") 
    results_writer.writerow(["RESULTS FOR INTERNAL STANDARDS..."])
    results_writer.writerow(["[sequence]", "[bio % backbone coverage]", "[bio % "+ion_type+" intensity]", "[Pearson r]", "[Fisher's z]", "[# of pairs]", "[z]"])
    
    #PCC: Calculate
    
    PCC_r_round=0
    PCC_list=[]
    query_PCCr = 0
    score_list=[]
    weighted_score_list=[]
    PCC_dud_ISP = []
    
    for i in range(0,len(top_bio_scans)):
        if len(top_bio_scans[i])>1 and len(top_syn_scans[i])>1:
            leading_bio_scan = deepcopy(top_bio_scans[i][3])
            leading_syn_scan = deepcopy(top_syn_scans[i][3])
            PCC_results = PCC_calculator(abund_thresh,PCC_abund_thresh,pro_mz_tol,leading_bio_scan,leading_syn_scan)
            PCC_r = PCC_results[1]
            biomin, synmin = PCC_results[3], PCC_results[4]
            if PCC_r != "Missing spectrum" and i == (len(top_bio_scans)-1):
                mirror_plot(hyphen, sequence, N_term_shift, C_term_shift, abund_thresh, biomin, synmin, top_bio_scans[len(top_bio_scans)-1][3], top_bio_scans[len(top_bio_scans)-1][4], top_bio_scans[len(top_bio_scans)-1][5], top_syn_scans[len(top_syn_scans)-1][3], top_syn_scans[len(top_syn_scans)-1][4], top_syn_scans[len(top_syn_scans)-1][5], ion_type, out_dir+"\\Figures\\", sequence_formatted + "_mirror.png")    
            if verbose == "Y" and PCC_r != "Missing spectrum" and i != (len(top_bio_scans)-1):
                ISP_sequence = top_bio_scans[i][0]
                ISP_sequence_formatted = ISP_sequence
                if N_term_shift_list[i] != 0:
                    ISP_sequence_formatted = "("+str(round(N_term_shift_list[i],2))+")"+ISP_sequence_formatted 
                if C_term_shift_list[i] != 0:
                    ISP_sequence_formatted = ISP_sequence_formatted+"("+str(round(C_term_shift_list[i],2))+")"
                mirror_plot(ISP_sequence.find("-"), ISP_sequence, N_term_shift_list[i], C_term_shift_list[i], abund_thresh, biomin, synmin, top_bio_scans[i][3], top_bio_scans[i][4], top_bio_scans[i][5], top_syn_scans[i][3], top_syn_scans[i][4], top_syn_scans[i][5], ion_type, out_dir+"\\Figures\\ISPs\\", ISP_sequence_formatted + "_mirror.png")                    
            if PCC_r != "Missing spectrum" and len(PCC_results[0]) < min_pairs_PCC:
                if i == (len(top_bio_scans)-1):
                    print("")
                    print("*** Too few pairs available for Pearson analysis ***")
                    results_writer.writerow([])
                    results_writer.writerow(["TOO FEW PAIRS AVAILABLE FOR PEARSON ANALYSIS"])
                    print()
                    results.close()
                    processing_time = timer(datetime.now(), start_time)
                    return(processing_time, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "TOO FEW PAIRS AVAILABLE FOR PEARSON ANALYSIS")
                else:
                    PCC_r = "Too few pairs available for Pearson analysis"
                    PCC_dud_ISP.append(top_bio_scans[i][0])
        elif len(top_bio_scans[i])<=1 and len(top_syn_scans[i])<=1:
            PCC_r = "Both the bio and syn run lack any spectra with correct precursor mass"
        elif len(top_bio_scans[i])<=1:
            PCC_r = "Bio run lacks any spectra with correct precursor mass"
        elif len(top_syn_scans[i])<=1:
            PCC_r = "Syn run lacks any spectra with correct precursor mass"
        if PCC_r == "Missing spectrum" and i == (len(top_bio_scans)-1):
            print("")
            print("*** In "+PCC_results[2]+", none of the spectra with precursor mass matching peptide of interest satisfy minimum scoring requirements ***")
            results_writer.writerow([])
            results_writer.writerow(["IN "+PCC_results[2].upper()+", NONE OF THE SPECTRA WITH PRECURSOR MASS MATCHING PEPTIDE OF INTEREST SATISFY MINIMUM SCORING REQUIREMENTS"])
            print()
            results.close()
            processing_time = timer(datetime.now(), start_time)
            return(processing_time, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "IN "+PCC_results[2].upper()+", NONE OF THE SPECTRA WITH PRECURSOR MASS MATCHING PEPTIDE OF INTEREST SATISFY MINIMUM SCORING REQUIREMENTS")
        if type(PCC_r)==float and i == (len(top_bio_scans)-1):
            query_PCCr =PCC_r
            syn_list_filtered = PCC_results[0] 
        elif type(PCC_r)==float and i != (len(top_bio_scans)-1):
            PCC_r_round=round(PCC_r,3)  
            pairs=len(PCC_results[0])  
            if PCC_r >= min_PCC:
                PCC_list.append(PCC_r)
            else:
                PCC_dud_ISP.append(top_bio_scans[i][0])
        else:
            PCC_r_round=PCC_r
        if i != (len(top_bio_scans)-1):
            if "lack" in str(PCC_r):
                print(f" {top_bio_scans[i][0]:24}  {PCC_r_round:69}")
                results_writer.writerow([top_bio_scans[i][0], PCC_r_round])
            elif PCC_r == "Missing spectrum":
                note = "In "+PCC_results[2]+", none of the spectra with correct precursor mass satisfy minimum scoring requirements"
                print(f" {top_bio_scans[i][0]:24}  {note:30}")
                results_writer.writerow([top_bio_scans[i][0], note])   
            elif PCC_r == "Too few pairs available for Pearson analysis":
                print(f" {top_bio_scans[i][0]:24}  {PCC_r_round:69}")
                results_writer.writerow([top_bio_scans[i][0], PCC_r_round])                
            else:
                Fisherz = atanh(PCC_r_round)
                print(f" {top_bio_scans[i][0]:24}  {round(top_bio_scans[i][1]):>23}  {round(top_bio_scans[i][2]):>21}  {PCC_r_round:>11}  {round(Fisherz,3):>12}  {pairs:>12}  {syn_charge_list[i]:>3}") 
                results_writer.writerow([top_bio_scans[i][0], round(top_bio_scans[i][1]), round(top_bio_scans[i][2]), PCC_r_round, round(Fisherz,3), pairs, syn_charge_list[i]])
                score_list.append(top_bio_scans[i][1])
                weighted_score_list.append(top_bio_scans[i][2])                    
    
    #Determine prediction interval for PCC distribution
    
    
    PCCr_threshold, query_PCCr_percentile = PCC_percentile_rank(PCC_list, percentile_thresh/100, query_PCCr, out_dir, sequence_formatted)
    
    if len(PCC_list) < min_intstd:
        PCC_outcome = "Too few internal standards"
    elif query_PCCr > PCCr_threshold and query_PCCr >= min_PCC:
        PCC_outcome = "PASS"
    else:
        PCC_outcome = "FAIL"
    
    query_PCC_syn, query_PCC_bio = [], []
    for i in range(0,len(syn_list_filtered)):
         query_PCC_syn.append(float(syn_list_filtered[i][1]))
         query_PCC_bio.append(float(syn_list_filtered[i][2]))
    
    #Generate plots
    
    regression_plot(query_PCC_syn, query_PCC_bio, "NA", "NA", "NA", "NA", sequence_formatted+" correlation", "syn raw intensity", "bio raw intensity", "Y", "Y", "PCC="+str(round(query_PCCr, 3)), out_dir+"\\Figures\\", sequence_formatted+"_correlation", "N", "N")
    swarm_plot(PCC_list, "standards", query_PCCr, sequence_formatted, "Pearson correlation coefficient (PCC)", "Pearson r", out_dir+"\\Figures\\", sequence_formatted, "_PCC_dist.png", PCCr_threshold, 1, round(query_PCCr_percentile,1))   
    
    #PCC: Results
    
    print()
    print("###############################################################################################################################",end="\n\n")
    print("PCC RESULTS FOR "+sequence+"...",end="\n\n")
    print("PCC:          ", round(query_PCCr, 3),    "      percentile:          ", str(round(query_PCCr_percentile,1))+"%", "      # of pairs:", len(syn_list_filtered))
    print("PCC threshold:", round(PCCr_threshold,3), "      percentile threshold:", str(round(percentile_thresh,1))+"%", end="\n\n")                 
    print("###############################################################################################################################",end="\n\n")
    
    results_writer.writerow([])
    results_writer.writerow(["PCC RESULTS FOR "+sequence+"..."])
    results_writer.writerow(["PCC:", round(query_PCCr, 3), "percentile:", str(round(query_PCCr_percentile,1))+"%", "# of pairs:", len(syn_list_filtered)])                   
    results_writer.writerow(["PCC threshold:", round(PCCr_threshold,3), "percentile threshold:", str(round(percentile_thresh,1))+"%"])
                             
    PCC_prep=open(str(out_dir + "\\Tables\\" + sequence_formatted + "_PCC_prep.csv"),"w",newline="")
    PCC_prep_writer=writer(PCC_prep)
    PCC_prep_writer.writerow(["Pearson_r", round(query_PCCr, 3), "peak_pairs", len(syn_list_filtered)])
    PCC_prep_writer.writerow([])
    PCC_prep_writer.writerow(["syn_m/z", "syn_abund", "bio_abund", "bio_m/z"])
    PCC_prep_writer.writerows(syn_list_filtered)
    PCC_prep.close()
       
    #Determine rough retention times (RTs) based on RT reported in MGF file for the best spectrum
    
    bio_roughRT = []
    bio_RTmzs = []
    for i in range(0, len(top_bio_scans)):
        if len(top_bio_scans[i]) == 1:
            bio_roughRT.append("Missing")
            bio_RTmzs.append("Missing")
        elif top_bio_scans[i][3] == "Missing spectrum" or (top_bio_scans[i][0] in PCC_dud_ISP):
            bio_roughRT.append("Missing")
            bio_RTmzs.append("Missing")
        else:
            RTINSECONDS = top_bio_scans[i][3][2][0]
            RT_sec = float(RTINSECONDS[12:])
            RT_min = RT_sec/60
            bio_roughRT.append(RT_min)
            bio_RTmzs.append((targets[i][2] + bio_charge_list[i]*1.0072764)/bio_charge_list[i])
         
    syn_roughRT = []
    syn_RTmzs = []
    for i in range(0, len(top_syn_scans)):
        if len(top_syn_scans[i]) == 1:
            syn_roughRT.append("Missing")
            syn_RTmzs.append("Missing")
        elif top_syn_scans[i][3] == "Missing spectrum":
            syn_roughRT.append("Missing")
            syn_RTmzs.append("Missing")
        else:
            RTINSECONDS = top_syn_scans[i][3][2][0]
            RT_sec = float(RTINSECONDS[12:])
            RT_min = RT_sec/60
            syn_roughRT.append(RT_min)
            syn_RTmzs.append((targets[i][2] + syn_charge_list[i]*1.0072764)/syn_charge_list[i])
    
    print("FILTERING MS1 FILE FOR BIOLOGICAL SAMPLE...", end="\n\n")    
    
    bio_MS1filtered = MS1_RTfilter(RTtol, bio_roughRT, bio_MS1path)
    
    print("FILTERING MS1 FILE FOR SYNTHETIC SAMPLE...", end="\n\n") 
    
    syn_MS1filtered = MS1_RTfilter(RTtol, syn_roughRT, syn_MS1path)
    
    print("DETERMINING RETENTION TIMES (IN MINUTES)...", end="\n\n")    
    print("[  sequence  ]            [bio RT]  [syn RT]  [syn RT - bio RT]", end="\n\n")    
    results_writer.writerow([])
    results_writer.writerow(["RETENTION TIMES (IN MINUTES)..."])
    results_writer.writerow(["[sequence]", "[bio RT]", "[syn RT]", "[syn RT - bio RT]"])
    
    bio_RTs, syn_RTs, delta_RTs = [], [], []
    for i in range(0, len(bio_RTmzs)):  
        if bio_RTmzs[i] != "Missing" and syn_RTmzs[i] != "Missing":
            bio_RT_results = extraction(bio_RTmzs[i], pre_mz_tol, bio_MS1filtered, bio_roughRT[i], RTtol)
            syn_RT_results = extraction(syn_RTmzs[i], pre_mz_tol, syn_MS1filtered, syn_roughRT[i], RTtol)
            bio_RTs.append(bio_RT_results[0])
            syn_RTs.append(syn_RT_results[0])
            delta_RTs.append(syn_RT_results[0] - bio_RT_results[0])  
            if i != len(bio_RTmzs) - 1:
                print(f" {top_bio_scans[i][0]:24}  {round(bio_RT_results[0],2):6}  {round(syn_RT_results[0],2):>8}  {round(syn_RT_results[0]-bio_RT_results[0],2):>17}") 
                results_writer.writerow([top_bio_scans[i][0], round(bio_RT_results[0],2), round(syn_RT_results[0],2), round(syn_RT_results[0]-bio_RT_results[0],2)])
                if verbose == "Y":
                    ISP_sequence = top_bio_scans[i][0]
                    ISP_sequence_formatted = ISP_sequence
                    if N_term_shift_list[i] != 0:
                        ISP_sequence_formatted = "("+str(round(N_term_shift_list[i],2))+")"+ISP_sequence_formatted 
                    if C_term_shift_list[i] != 0:
                        ISP_sequence_formatted = ISP_sequence_formatted+"("+str(round(C_term_shift_list[i],2))+")"
                    EIC(bio_RT_results[1], bio_RT_results[2], bio_RT_results[0], "biological run", syn_RT_results[1], syn_RT_results[2], syn_RT_results[0], "synthetic run", ISP_sequence_formatted, out_dir+"\\Figures\\ISPs\\", ISP_sequence_formatted + "_EIC.png")

            else:
                bio_EICx = bio_RT_results[1]
                bio_EICy = bio_RT_results[2]
                syn_EICx = syn_RT_results[1]
                syn_EICy = syn_RT_results[2]
    
    EIC(bio_EICx, bio_EICy, bio_RTs[len(bio_RTs)-1], "biological run", syn_EICx, syn_EICy, syn_RTs[len(syn_RTs)-1], "synthetic run", sequence_formatted, out_dir+"\\Figures\\", sequence_formatted + "_EIC.png")
    
    regression_plot(bio_RTs[0:len(bio_RTs)-1], syn_RTs[0:len(bio_RTs)-1], "standards", bio_RTs[len(bio_RTs)-1], syn_RTs[len(bio_RTs)-1], sequence_formatted, "retention time (RT) in minutes", "RT in biological run", "RT in synthetic run", "N", "N", "N", out_dir+"\\Figures\\", sequence_formatted+"_RT", min_RT, max_RT)
    regression_plot(bio_RTs[0:len(bio_RTs)-1], delta_RTs[0:len(bio_RTs)-1], "standards", bio_RTs[len(bio_RTs)-1], delta_RTs[len(bio_RTs)-1], sequence_formatted, "raw delta RT", "RT in biological run (minutes)", "syn RT - bio RT (minutes)", "N", "N", "N", out_dir+"\\Figures\\", sequence_formatted+"_deltaRT", min_RT, max_RT)
    
    ref_RTs, test_RTs = [],[]
    for i in range(0, len(bio_RTs)-1):
        if min_RT <= bio_RTs[i] <= max_RT:
            ref_RTs.append(bio_RTs[i])
            test_RTs.append(syn_RTs[i])
    
    query_ref_RT, query_test_RT = bio_RTs[len(bio_RTs)-1], syn_RTs[len(syn_RTs)-1]
    
    if query_ref_RT >= max(ref_RTs) or query_ref_RT <= min(ref_RTs):
        RT_outcome = "RT for peptide of interest in the biological sample is outside the range of the internal standards"
        print()
        print("###############################################################################################################################",end="\n\n")
        print("RT FOR PEPTIDE OF INTEREST IN THE BIOLOGICAL SAMPLE IS OUTSIDE THE RANGE OF THE INTERNAL STANDARDS",end="\n\n")
        print("###############################################################################################################################",end="\n\n")
        results_writer.writerow([])
        results_writer.writerow(["RT FOR PEPTIDE OF INTEREST IN THE BIOLOGICAL SAMPLE IS OUTSIDE THE RANGE OF THE INTERNAL STANDARDS"])
    else:
        sorted_RTs, RT_pred_deltas, query_RT_pred_delta = RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT)
        
        delta_lo, delta_hi, query_percentile = RT_percentile_rank(RT_pred_deltas, percentile_thresh/100, query_RT_pred_delta, out_dir, sequence_formatted)
        
        if delta_lo > -manual_RTdev_thresh:
            delta_lo = -manual_RTdev_thresh
        if delta_hi < manual_RTdev_thresh:
            delta_hi = manual_RTdev_thresh

        if len(RT_pred_deltas) < min_intstd:
            RT_outcome = "Too few internal standards"
        elif delta_hi > query_RT_pred_delta > delta_lo:
            RT_outcome = "PASS"
        else:
            RT_outcome = "FAIL"
        
        print()
        print("###############################################################################################################################",end="\n\n")
        print("RT RESULTS (IN MINUTES) FOR "+sequence+"...",end="\n\n")
        print("raw delta RT:", round(syn_RT_results[0]-bio_RT_results[0],3), "      deviation from expected RT:", round(query_RT_pred_delta,3),                           "               percentile:          ", str(round(query_percentile, 1))+"%")
        print(                                           "                          acceptable deviation:      ", str(round(delta_lo, 3))+" to "+str(round(delta_hi, 3)),          "      percentile threshold:", str(round(percentile_thresh,1))+"%", end="\n\n")                    
        print("###############################################################################################################################",end="\n\n")
        results_writer.writerow([])
        results_writer.writerow(["RT RESULTS (IN MINUTES) FOR "+sequence+"..."])
        results_writer.writerow(["raw delta RT:", round(syn_RT_results[0]-bio_RT_results[0],3), "deviation from expected RT:", round(query_RT_pred_delta,3), "percentile:", str(round(query_percentile, 1))+"%"])
        results_writer.writerow(["", "", "acceptable deviation:", str(round(delta_lo, 3))+" to "+str(round(delta_hi, 3)), "percentile threshold:", str(round(percentile_thresh,1))+"%"])
        
#        sorted_RTs_bio = []
#        for i in range(1, len(sorted_RTs)-1):
#            sorted_RTs_bio.append(sorted_RTs[i][0])
#        regression_plot(sorted_RTs_bio, RT_pred_deltas, "standards", query_ref_RT, query_RT_pred_delta, sequence_formatted, "deviation from expected RT", "RT in biological run (minutes)", "expected RT - syn RT (minutes)", "N", "N", "N", out_dir+"\\Figures", sequence_formatted+"_expectedRT")
        
        swarm_plot(RT_pred_deltas, "standards", query_RT_pred_delta, sequence_formatted, "deviation from expected RT", "expected RT - syn RT (minutes)", out_dir+"\\Figures\\", sequence_formatted, "_expectedRT_dist.png", delta_lo, delta_hi, round(query_percentile,1))   
    
    print("PCC test: "+PCC_outcome+"     RT test: "+RT_outcome, end="\n\n")
    processing_time = timer(datetime.now(), start_time)
    print("Elapsed time (min:sec):",processing_time,end="\n\n")
    results_writer.writerow([])
    results_writer.writerow(["PCC test:", PCC_outcome, "RT test:", RT_outcome])
    results_writer.writerow([])
    results_writer.writerow(["processing time (min:sec):", processing_time])
    results.close()
    print("Your files are ready in the folder " + out_dir,end="\n\n")
    print("###############################################################################################################################",end="\n\n")

    if RT_outcome == "RT for peptide of interest in the biological sample is outside the range of the internal standards":
        return(processing_time, round(query_PCCr, 3), str(round(query_PCCr_percentile,1))+"%", "N/A", "N/A", PCC_outcome, "N/A", RT_outcome)
    else:
        return(processing_time, round(query_PCCr, 3), str(round(query_PCCr_percentile,1))+"%", round(query_RT_pred_delta,2), str(round(query_percentile, 1))+"%", PCC_outcome, RT_outcome)

##################################################################################################################################################################
