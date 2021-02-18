'''
is it merging? if so, write your own merging method in which you sum intensities and average the m/zs; could use PCC_calculator as basis for this
you would also basically keep a running list of all precursors with m/z and RT. 
with each new one, you look back at the list and once you are in the right RT range see if any of the previous scans have the same m/z
if they do, could then check for some level of similarity and if they pass that then you merge the scans
NOTE: the analyses seem to work just fine on the benchmarking data, so what is different about that data and the hip data that doesn't work? i think it's the abundance of the bio hip maybe
i can still use spectrum mill mzxml output as a starting point with the version of program that has msconvert front end, just have .d and .mzxml in folder and specify different file types in msconvert_settings

######have program spit out the ms1 list it makes that it uses for precursor refinement so you can see if it looks right

need to add a note about keeping low abundance peaks for precursor refine in method settings

change log:
    pro_mz_tol doubled for PCC
    version reported in detailed summary
    wrote my own precursor refine
    
other big changes i want to make: 
    removing some of the redundancy in settings, like getting rid of abundance thresh and just letting msconvert take care of that
    making msconvert setttings kind of fixed so don't have to understand how to enter CLI arguments and merge msconvert_settings and settings files into one file
'''

from os import rename
from csv import reader, writer
from copy import deepcopy

##################################################################################################################################################################

def precursor_refine(sample_ms1, sample_mgf):
        
    MS1=open(sample_ms1)
    MS1content=reader(MS1, delimiter=" ")
    MS1_data = [[],[]]
    scan_start = 5
    scan = []
    for row in MS1content:
        row = list(row)
        if "S\t0\t0" in row:
            if MS1content.line_num > 5:
                MS1_data[1].append(scan)
            scan_start = int(MS1content.line_num) 
            scan = []  
        elif "RTime" in row[0]:
            RT_minutes = float(row[0][8:]) 
            MS1_data[0].append(RT_minutes)
        elif (MS1content.line_num > (scan_start + 4)):
                scan.append(float(row[0]))
    MS1.close()

    unrefined_file_name = sample_mgf[0:len(sample_mgf)-4] + "_raw.mgf"
    rename(sample_mgf, unrefined_file_name)    
    
    raw = open(unrefined_file_name)
    contents = reader(raw, delimiter=" ")

    refined_file_name = sample_mgf[0:len(sample_mgf)-4] + "_precursors_refined.mgf"
    refined = open(refined_file_name,"w",newline="")
    refined_writer = writer(refined,delimiter=" ")

    for row in contents:
        info = list(row)
        if "RTINSECONDS" in info[0]:
            RTINSECONDS=str(info[0])
            RT=float(RTINSECONDS[12:])/60
        elif "PEPMASS" in info[0]:
            PEPMASS = str(info[0])
            obs_pre_mz = float(PEPMASS[8:])
            if RT > MS1_data[0][0]:     #in case ms1 scan preceding the ms2 got removed by msconvert
                for i in range(0, len(MS1_data[0])):
                    if MS1_data[0][i] > RT:
                        for j in range(0, len(MS1_data[1][i-1])):
                            if MS1_data[1][i-1][j] > obs_pre_mz:
                                left = MS1_data[1][i-1][j-1]
                                right = MS1_data[1][i-1][j]
                                left_dif = obs_pre_mz - left
                                right_dif = right - obs_pre_mz
                                if left_dif < right_dif:
                                    refined_pre_mz = left
                                elif right_dif < left_dif:
                                    refined_pre_mz = right
                                info[0] = "PEPMASS=" + str(refined_pre_mz)
                                break
                        break     
        refined_writer.writerow(info)

    raw.close()
    refined.close()

##################################################################################################################################################################
    
def spectrum_merge(refined_mgf, RTtol, pre_mz_tol, pro_mz_tol):
    
    pre_mz_tol = 2*pre_mz_tol
    pro_mz_tol = 2*pro_mz_tol
    
    def merge(reference_scan, scan, pro_mz_tol): #Keeps RT from first scan, averages pre_mz; for product ions, adds together using pro_mz_tol to determine if peaks are the same
        
        merged_spectrum = []
        merged_spectrum.append(reference_scan[0])
        merged_spectrum.append(reference_scan[1])  
        merged_spectrum.append(reference_scan[2]) 

        reference_pre_mz = float(reference_scan[3][0][8:])
        pre_mz = float(scan[3][0][8:])
        merged_mz = (reference_pre_mz + pre_mz)/2
        merged_mz_statement = "PEPMASS=" + str(merged_mz)
        merged_abundance = float(reference_scan[3][1]) + float(scan[3][1]) 
        merged_spectrum.append([merged_mz_statement, merged_abundance])

        merged_spectrum.append(reference_scan[4]) 
        ref_peaks = []
        for i in range(5, len(reference_scan) - 1):
            peak = [float(reference_scan[i][0]), float(reference_scan[i][1])]
            ref_peaks.append(peak)        
        peaks = []
        for i in range(5, len(scan) - 1):
            peak = [float(scan[i][0]), float(scan[i][1])]            
            peaks.append(peak)
        
        for i in range(0, len(peaks)):
            for j in range(0, len(ref_peaks)):
                if ((ref_peaks[j][0]-ref_peaks[j][0]/1000000*pro_mz_tol) <= peaks[i][0] <= (ref_peaks[j][0]+ref_peaks[j][0]/1000000*pro_mz_tol)):
                    ref_peaks[j][0] = (ref_peaks[j][0] + peaks[i][0])/2
                    ref_peaks[j][1] = ref_peaks[j][1] + peaks[i][1]
                    break
                elif ((peaks[i][0] - ref_peaks[j][0])/ref_peaks[j][0]*1000000) > pro_mz_tol:
                    ref_peaks.append(peaks[i])
                    break        

        ref_peaks.sort()   
        for i in range(0, len(ref_peaks)):
            merged_spectrum.append(ref_peaks[i])
        
        merged_spectrum.append(["END","IONS"])
        
        return merged_spectrum
    
    unmerged = open(refined_mgf)
    contents = reader(unmerged, delimiter=" ")

    merged_file_name = refined_mgf[0:len(refined_mgf)-23] + ".mgf"
    merged = open(merged_file_name,"w",newline="")
    merged_writer = writer(merged,delimiter=" ")

    MS2_scans = []
    for row in contents:
        if "BEGIN" in row:
            scan_start=int(contents.line_num)
            scan=[["BEGIN","IONS"]]
        elif contents.line_num == (scan_start + 1):
            scan.append(list(row))
        elif contents.line_num == (scan_start + 2):
            RTINSECONDS=list(row)
            scan.append(RTINSECONDS)
            RTINSECONDS=str(RTINSECONDS[0])
            RT=float(RTINSECONDS[12:])/60
        elif contents.line_num == (scan_start + 3): 
            PEPMASS=list(row)
            scan.append(PEPMASS)
            PEPMASS=PEPMASS[0]
            obs_pre_mz=float(PEPMASS[8:])
        elif contents.line_num == (scan_start + 4):
            CHARGE=list(row)
            scan.append(CHARGE)
            CHARGE=str(CHARGE[0])
            obs_z=int(CHARGE[7])        
        elif (contents.line_num > (scan_start + 4)) and ("END" not in row):
            peak=list(row)
            scan.append(peak)
        elif "END" in row:
            scan.append(list(row))
            if MS2_scans == []:
                MS2_scans.append(scan)
            else:
                for i in range(len(MS2_scans)-1, -1, -1):
                    reference_RT = float(MS2_scans[i][2][0][12:])/60
                    reference_pre_mz = float(MS2_scans[i][3][0][8:])
                    reference_z = int(MS2_scans[i][4][0][7])
                    if RT - reference_RT > RTtol:
                        MS2_scans.append(scan)
                        break
                    elif (RT - reference_RT <= RTtol) and (reference_z == obs_z) and ((reference_pre_mz-reference_pre_mz/1000000*pre_mz_tol) <= obs_pre_mz <= (reference_pre_mz+reference_pre_mz/1000000*pre_mz_tol)): 
                        reference_scan = deepcopy(MS2_scans[i])
                        merged_spectrum = merge(reference_scan, scan, pro_mz_tol)
                        MS2_scans[i] = merged_spectrum
                        break
                                           
    for i in range(0, len(MS2_scans)):
        merged_writer.writerows(MS2_scans[i])

    unmerged.close()
    merged.close()
    
##################################################################################################################################################################
