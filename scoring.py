import copy

##################################################################################################################################################################

def LR_score(abund_thresh,pro_mz_tol,min_score,min_weighted_score,hits,hits_writer,top_hit_writer, fixed_charge, ion_type):

    ##Put explanation here of the data structure (really the return value) (so show the nested lists and what is in each position) and detailed explanations of how the program is running step by step
##    [
##    ["YSAHE-EHHYDK", leading_score, leading_weighted_score, leading_scan, leading_L_ions, leading_R_ions],
##    ["HEHISS-DYAGK", leading_score, leading_weighted_score, leading_scan, leading_L_ions, leading_R_ions],
##    ...
##    ]

    query_details = []
    top_hits=[]
    charge_list = []
    for i in range(0,len(hits)):
        hyphen = hits[i][1].find("-")
        if hits[i][0] == "putative":
            print()
            print("For each spectrum with the correct precursor mass:")
            print()
            if hyphen != -1:
                print("[L pep backbone coverage]  [R pep backbone coverage]  [% backbone coverage]  [% "+ion_type+" intensity]  [z]")
                print()
            else:
                print("[% backbone coverage]  [% "+ion_type+" intensity]  [z]")
                print()
        top_hits.append([hits[i][1]])
        leading_charge = "N/A"
        if len(hits[i])>8:
            L_bonds=hits[i][4]
            L_ions=hits[i][5]
            R_bonds=hits[i][6]
            R_ions=hits[i][7]
            leading_score, leading_weighted_score=0,0
            chosen_L_score, chosen_R_score = "N/A", "N/A"
            scans_passed=0
            leading_L_ions, leading_R_ions=[],[] 
            leading_scan="Missing spectrum"
            for j in range(8,len(hits[i])):
                L_total, R_total=0,0
                L_score, R_score, score="N/A", "N/A", 0
                total_signal, matched_signal, weighted_score = 0,0,0
                junction_b, junction_y = "no", "no"
                RTINSECONDS=hits[i][j][2]
                RTINSECONDS=str(RTINSECONDS[0])
                RT=float(RTINSECONDS[12:])/60
                PEPMASS = hits[i][j][3]
                PEPMASS = PEPMASS[0]
                obs_pre_mz = float(PEPMASS[8:])
                upper_premz, lower_premz = (obs_pre_mz + obs_pre_mz/1000000*pro_mz_tol), (obs_pre_mz - obs_pre_mz/1000000*pro_mz_tol)
                CHARGE=hits[i][j][4]
                CHARGE=str(CHARGE[0])
                obs_z=int(CHARGE[7])   
                for k in range(5,len(hits[i][j])-1):
                    peak=(hits[i][j][k])
                    ###You aren't accounting for the possibility that a peak could match more than one predicted ion and more than one peak could match a predicted ion. Could use PCCprep approach.
                    if (float(peak[1])) > abund_thresh and (float(peak[0]) > upper_premz or float(peak[0]) < lower_premz):
                        total_signal=total_signal+float(peak[1])
                        for l in range (0, len(L_ions)):
                            for m in range(1,5):
                                if ((L_ions[l][m]-L_ions[l][m]/1000000*pro_mz_tol) <= float(peak[0]) <= (L_ions[l][m]+L_ions[l][m]/1000000*pro_mz_tol)):
                                    L_ions[l][5]=L_ions[l][5]+1
                                    L_ions[l][m+6]=peak[0]
                                    L_ions[l][m+11]=peak[1]
                        if hyphen != -1:
                            for l in range (0, len(R_ions)):
                                for m in range(1,5):
                                    if ((R_ions[l][m]-R_ions[l][m]/1000000*pro_mz_tol) <= float(peak[0]) <= (R_ions[l][m]+R_ions[l][m]/1000000*pro_mz_tol)):
                                        R_ions[l][5]=R_ions[l][5]+1
                                        R_ions[l][m+6]=peak[0]
                                        R_ions[l][m+11]=peak[1]
                for n in range(len(L_ions)):
                    if L_ions[n][5]>0:
                        L_total=L_total+1
                        if n == (len(L_ions) - 1):
                            junction_b = "yes"
                    for o in range(12,16):
                        if L_ions[n][o] != "":
                            matched_signal = matched_signal + float(L_ions[n][o])   
                if hyphen != -1:
                    for n in range(len(R_ions)):
                        if R_ions[n][5]>0:
                            R_total=R_total+1
                            if n == 0:
                                junction_y = "yes"
                        for o in range(12,16):
                            if R_ions[n][o] != "":
                                matched_signal = matched_signal + float(R_ions[n][o])
                    if junction_b == "yes" and junction_y == "yes":
                        score = (L_total + R_total - 1)/(L_bonds + R_bonds -1)*100 #addresses issue of duplicated junction bond for hybrids                
                    else: 
                        score = (L_total + R_total)/(L_bonds + R_bonds -1)*100
                else:
                    score = (L_total)/(L_bonds)*100 
                if total_signal > 0:
                    weighted_score = matched_signal/total_signal*100
                else:
                    weighted_score = 0

                ##hits_writer.writerow(["SCAN #",scan_num])
                if hits[i][0] == "putative":    
                    hits_writer.writerow(["RT (minutes)", RT, "charge", obs_z])                    
                    hits_writer.writerow([])
                    hits_writer.writerow(["pred_L_ions","","","","","","","obs_m/z","","","","","obs_abund"])
                    if ion_type == "b/y":
                        hits_writer.writerow(["bond #", "b+", "b++", "y+", "y++", "matches", "", "b+", "b++", "y+", "y++", "", "b+", "b++", "y+", "y++"])
                    elif ion_type == "c/z":
                        hits_writer.writerow(["bond #", "c+", "c++", "z+", "z++", "matches", "", "c+", "c++", "z+", "z++", "", "c+", "c++", "z+", "z++"])                        
                    hits_writer.writerows(L_ions)
                    hits_writer.writerow([])
                    if hyphen != -1:
                        hits_writer.writerow(["pred_R_ions","","","","","","","obs_m/z","","","","","obs_abund"])
                        if ion_type == "b/y":
                            hits_writer.writerow(["bond #", "b+", "b++", "y+", "y++", "matches", "", "b+", "b++", "y+", "y++", "", "b+", "b++", "y+", "y++"])
                        elif ion_type == "c/z":
                            hits_writer.writerow(["bond #", "c+", "c++", "z+", "z++", "matches", "", "c+", "c++", "z+", "z++", "", "c+", "c++", "z+", "z++"])                                                   
                        hits_writer.writerows(R_ions)
                        hits_writer.writerow([])
                        L_score = str(L_total)+"/"+str(L_bonds)
                        R_score = str(R_total)+"/"+str(R_bonds)
                        hits_writer.writerow(["L pep backbone coverage", "'"+L_score, "R pep backbone coverage", "'"+R_score, "% backbone coverage =", round(score,2), "% "+ion_type+" intensity =", round(weighted_score,2)])
                        print(f"{L_score:>24}   {R_score:>24}  {round(score):>21}  {round(weighted_score):>17}  {obs_z:>3}")
                        query_details.append([L_score, R_score, score, weighted_score, obs_z])
                    else:
                        hits_writer.writerow(["% backbone coverage =", round(score,2), "% "+ion_type+" intensity =", round(weighted_score,2)]) 
                        print(f"{round(score):>20}  {round(weighted_score):>17}  {obs_z:>3}")
                        query_details.append([score, weighted_score, obs_z])
                    hits_writer.writerow([])
                    hits_writer.writerow([])
                    hits_writer.writerow(["********","********","********","********","********","********","********","********","********","********","********","********","********","********","********","********"])
                    hits_writer.writerow([])
                    hits_writer.writerow([])
                    scans_passed=scans_passed+1
                if fixed_charge == "N" or obs_z == fixed_charge[i] or fixed_charge[i] == "N/A":
                    if score>leading_score and score >= min_score and weighted_score >= min_weighted_score:
                        leading_score=score
                        leading_weighted_score=weighted_score
                        chosen_L_score, chosen_R_score = L_score, R_score
                        leading_scan=copy.deepcopy(hits[i][j])
                        leading_L_ions=copy.deepcopy(L_ions) 
                        leading_R_ions=copy.deepcopy(R_ions)
                        leading_charge = obs_z
                    elif score==leading_score and score >= min_score and weighted_score >= min_weighted_score:
                        if weighted_score > leading_weighted_score:
                            leading_score=score
                            leading_weighted_score=weighted_score
                            chosen_L_score, chosen_R_score = L_score, R_score
                            leading_scan=copy.deepcopy(hits[i][j])
                            leading_L_ions=copy.deepcopy(L_ions) 
                            leading_R_ions=copy.deepcopy(R_ions)
                            leading_charge = obs_z
                for p in range (len(L_ions)):
                    L_ions[p][5:]=[0,"","","","","","","","","",""]
                for p in range (len(R_ions)):
                    R_ions[p][5:]=[0,"","","","","","","","","",""]
            top_hits[i].append(leading_score)
            top_hits[i].append(leading_weighted_score)
            top_hits[i].append(leading_scan)
            top_hits[i].append(leading_L_ions)
            top_hits[i].append(leading_R_ions)
        charge_list.append(leading_charge)

#do i want to have winning spectrum chosen by weighted score?
            
    print()
    if fixed_charge != "N":
        print("NOTE: Algorithm chooses winning spectrum only from those spectra with a corresponding precursor charge that matches the")
        print("      precursor charge for the winning spectrum chosen from the biological sample.")
        print()
    print("scans with correct precursor mass:",scans_passed,"     winning % backbone coverage:",round(leading_score,2), "      winning % "+ion_type+" intensity:",round(leading_weighted_score,2), "     z:",leading_charge, end="\n\n")
    top_hit_writer.writerows(leading_scan)
    
    return(top_hits, charge_list, chosen_L_score, chosen_R_score, query_details)
##Note you changed from leading scan to top hits as the return, no need to translate this into the pcc calculator and have pcc calculator cycle through all top hits for syn and bio, but also have to do this with the filter
##first as well

##################################################################################################################################################################