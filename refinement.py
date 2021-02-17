'''
run the new data (and can try with old too) using msconvert-converted vs specmill-converted; how do the results differ?
if behaves like you expect, figure out what specmill does that helps
is it precursor refine? if so, see how differently precursor refine works for the two programs; do i need to write my own precursor refine? i think you would just pick the peak closest in mass and that is biggest
why doesn't vendor software wait and do the mass accuracy calibration after everything is done so it could easily know which parent peak the ms2 metadata precursor corresponds to and just change it with that?
It's because the quad tried to grab that spot, so you basically just want to pick whatever parent ion falls closest to where the quad targeted!!!!!!!!!!!!!!!!
is it merging? if so, write your own merging method in which you sum intensities and average the m/zs; could use PCC_calculator as basis for this
you would also basically keep a running list of all precursors with m/z and RT. 
with each new one, you look back at the list and once you are in the right RT range see if any of the previous scans have the same m/z
if they do, could then check for some level of similarity and if they pass that then you merge the scans
NOTE: the analyses seem to work just fine on the benchmarking data, so what is different about that data and the hip data that doesn't work? i think it's the abundance of the bio hip maybe
i can still use spectrum mill mzxml output as a starting point with the version of program that has msconvert front end, just have .d and .mzxml in folder and specify different file types in msconvert_settings




change log:
    pro_mz_tol doubled for PCC
    version reported in detailed summary
    wrote my own precursor refine
    
other big changes i want to make: 
    removing some of the redundancy in settings, like getting rid of abundance thresh and just letting msconvert take care of that
    making msconvert setttings kind of fixed so don't have to understand how to enter CLI arguments and merge msconvert_settings and settings files into one file
    need to close some file or something in script so that when running in spyder it doens't keeping me from moving files saying that is in use by other program
'''

from os import rename
from csv import reader, writer

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

    refined = open(sample_mgf,"w",newline="")
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
    
def spectrum_merge():
    return()
    
##################################################################################################################################################################
