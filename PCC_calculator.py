from math import sqrt

##################################################################################################################################################################

def PCC_calculator(abund_thresh,PCC_abund_thresh,pro_mz_tol,leading_bio_scan,leading_syn_scan):

    '''    
    PCC: Pair data
        For each peak in the syn spectrum, check the peaks in the bio spectrum one by one to see if any have an m/z within +/-pro_mz_tol. If so, pair the two peaks.
        If more than one match is found, use the one with the closest m/z.
        If a match cannot be found, designate an abundance of zero for this m/z in the bio spectrum.
    PCC: Filter noise from peak lists (based on fixed, user-defined threshold. Could use alternative approaches like avg + 2SD)
    '''

    warning = "none"
    syn_list_filtered = []
    biomin, synmin = 0, 0
    if leading_bio_scan==("Missing spectrum") or leading_syn_scan==("Missing spectrum"):
        PCC_r="Missing spectrum"
        if leading_bio_scan==("Missing spectrum") and leading_syn_scan==("Missing spectrum"):
            warning = "both runs"
        elif leading_bio_scan==("Missing spectrum"):
            warning = "bio run"
        else:
            warning = "syn run"
        
    else:
        bio_pre = float(leading_bio_scan[3][0][8:])
        z = int(leading_bio_scan[4][0][7])
        syn_pre = float(leading_syn_scan[3][0][8:])
        increment = 1/z
        bio_pre_list = [bio_pre, bio_pre + increment, bio_pre + 2*increment]
        syn_pre_list = [syn_pre, syn_pre + increment, syn_pre + 2*increment]
        bio_window, syn_window = [], []
        for i in range(0, len(bio_pre_list)):
            bio_window.append([bio_pre_list[i] - bio_pre_list[i]/1000000*pro_mz_tol, bio_pre_list[i] + bio_pre_list[i]/1000000*pro_mz_tol])
        for i in range(0, len(syn_pre_list)):
            syn_window.append([syn_pre_list[i] - syn_pre_list[i]/1000000*pro_mz_tol, syn_pre_list[i] + syn_pre_list[i]/1000000*pro_mz_tol])
        
        bio_list=[]
        for i in range(5,len(leading_bio_scan)-1):
            bio_list.append(leading_bio_scan[i])
        syn_list=[]
        for i in range(5,len(leading_syn_scan)-1):
            syn_list.append(leading_syn_scan[i])
        for i in range(0,len(syn_list)):
            mz=float(syn_list[i][0])
            for j in range(0,len(bio_list)):
                if ((mz-mz/1000000*pro_mz_tol) <= float(bio_list[j][0]) <= (mz+mz/1000000*pro_mz_tol)) and (len(syn_list[i])==2):
                    syn_list[i].append(bio_list[j][1])
                    syn_list[i].append(bio_list[j][0])
                elif ((mz-mz/1000000*pro_mz_tol) <= float(bio_list[j][0]) <= (mz+mz/1000000*pro_mz_tol)) and (len(syn_list[i])>2) and (abs(mz-float(bio_list[j][0])) < abs(mz-float(syn_list[i][3]))):
            #If the previous match and the current match are equally good, algorithm keeps the first match (lower m/z).
            #This makes sense because software is scanning from low to high so you would prioritize what is likely to be the C12 peak.
                    syn_list[i][2]=bio_list[j][1]
                    syn_list[i][3]=bio_list[j][0]
                elif len(syn_list[i])==2 and j==(len(bio_list)-1):
                    syn_list[i].append(0)
                    syn_list[i].append(0)
            #Find any instances where a bio value was used more than once.
        for k in range(0,len(bio_list)):
            times_used=0
            uses=[]
            mz=float(bio_list[k][0])
            for l in range(0,len(syn_list)):
                if float(syn_list[l][3])==mz:
                    uses.append(syn_list[l][0])
                    times_used=times_used+1
            if times_used>1:
                uses_diff=[]
                for m in range(0,len(uses)):
                    uses_diff.append(abs(float(uses[m])-mz))
                uses_diff_sort=[]
                uses_diff_sort=uses_diff_sort+uses_diff
                uses_diff_sort.sort()
                best=uses_diff_sort[0]
                best_mz=uses[uses_diff.index(best)]
            #If more than one match had the same diff, index reports the one with the lowest m/z, and that is the one that will be used.
                for n in range(0,len(syn_list)):
                    if (float(syn_list[n][3])==mz) and (float(syn_list[n][0])!=float(best_mz)):
                        syn_list[n][2]=0
                        syn_list[n][3]=0
            #Now that you deleted this match, you need to see if there was a rank2 match in the bio spectrum. Make sure you don't add back the rank1 match.
                        syn_mz=float(syn_list[n][0])
                        for o in range(0,len(bio_list)):
                            if ((syn_mz-syn_mz/1000000*pro_mz_tol) <= float(bio_list[o][0]) <= (syn_mz+syn_mz/1000000*pro_mz_tol)) and (float(bio_list[o][0]) != mz) and (abs(syn_mz-float(bio_list[o][0])) < abs(syn_mz-float(syn_list[n][3]))):
            #The syn_list row here will always be 4 columns long so don't need to check its length. 
                                syn_list[n][2]=bio_list[o][1]
                                syn_list[n][3]=bio_list[o][0]
            #During this process, you might have made a second instance of an bio m/z you already checked, so you need to report an error if this happens...
            #I realized that whenever i have it put in the rank2 match to replace a rank 1 that had already been used, if I just have the criteria that the rank 2 can't already exist as a hit then that would take care of it.
            #I don't know though if that approach would ensure that this "rank2" match ended up assigned to the appropriate syn m/z, so maybe i have to do a nested check here where if it does already exist elsewhere it decides which spot is better                        
        for p in range(0,len(bio_list)):
            times_used=0
            uses=[]
            mz=float(bio_list[p][0])
            for q in range(0,len(syn_list)):
                if float(syn_list[q][3])==mz:
                    uses.append(syn_list[q][0])
                    times_used=times_used+1
            if times_used>1:
                print("")
                print("The following m/z in the bio spectrum was used more than once:")
                print(mz)
            elif times_used==0:
                syn_list.append([0,0,bio_list[p][1],bio_list[p][0]])
            #Check if you have any instances where when matching bio peaks to syn peaks the m/z got out of order.
        for r in range(0,len(syn_list)-1):
            if (float(syn_list[r][3]) > float(syn_list[r+1][3])) and (float(syn_list[r+1][3]) != 0) and (float(syn_list[r+1][0]) != 0):
                print("")
                print("The following peaks from the bio file were switched during matching:")
                print(syn_list[r][3])
                print(syn_list[r+1][3])

        #PCC: Filter
        
        #convert pcc abund thresh into a number relative to max and then add to conditions below
        
        syn_list_nopre = []
        for i in range(len(syn_list)):            
            precursor = "no"
            for j in range(0, len(bio_window)):
                if bio_window[j][0] <= float(syn_list[i][3]) <= bio_window[j][1]:
                    precursor = "yes"
            for j in range(0, len(syn_window)):
                if syn_window[j][0] <= float(syn_list[i][0]) <= syn_window[j][1]:
                    precursor = "yes"
            if precursor == "no":
                syn_list_nopre.append(syn_list[i])

        syn1max = max(syn_list_nopre, key = lambda x: float(x[1]))
        syn2max = max(syn_list_nopre, key = lambda x: float(x[2]))
        synmin = float(syn1max[1])*PCC_abund_thresh/100
        biomin = float(syn2max[2])*PCC_abund_thresh/100

        for i in range(len(syn_list_nopre)):            
            if (float(syn_list_nopre[i][1])>abund_thresh and float(syn_list_nopre[i][1])>synmin) or (float(syn_list_nopre[i][2])>abund_thresh and float(syn_list_nopre[i][2])>biomin):
                syn_list_filtered.append(syn_list_nopre[i])

        #PCC: Calculate sample PCC

        filtered_bio_total, filtered_syn_total = 0, 0
        for i in range(0, len(syn_list_filtered)):
            bio_abund=float(syn_list_filtered[i][2])
            syn_abund=float(syn_list_filtered[i][1])
            filtered_bio_total=filtered_bio_total+bio_abund
            filtered_syn_total=filtered_syn_total+syn_abund
        filtered_bio_avg=filtered_bio_total/len(syn_list_filtered)
        filtered_syn_avg=filtered_syn_total/len(syn_list_filtered)
        r_numerator=0
        for i in range(0, len(syn_list_filtered)):
            r_numerator=r_numerator+(float(syn_list_filtered[i][2])-filtered_bio_avg)*(float(syn_list_filtered[i][1])-filtered_syn_avg)
        sum_bio_sq=0
        sum_syn_sq=0
        for i in range(0, len(syn_list_filtered)):
            sum_bio_sq=sum_bio_sq+(float(syn_list_filtered[i][2])-filtered_bio_avg)**2
            sum_syn_sq=sum_syn_sq+(float(syn_list_filtered[i][1])-filtered_syn_avg)**2
        PCC_r=r_numerator/(sqrt(sum_bio_sq)*sqrt(sum_syn_sq))

    return(syn_list_filtered,PCC_r,warning,biomin,synmin)
    
##################################################################################################################################################################
