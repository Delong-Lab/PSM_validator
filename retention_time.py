from csv import reader

##################################################################################################################################################################

'''
Earlier in analysis, the highest-scoring spectrum in each MGF file was chosen for each peptide sequence.
The RT for the spectrum was recorded as the rough RT.
This function filters the MS1 files to contain only those spectra within a certain time window (RTtol) around each RT. 
This greatly reduces downstream processing time.
'''
    
def MS1_RTfilter(RTtol, roughRT_list, MS1path):

    RTbounds = [[],[]]
    for i in range(0, len(roughRT_list)):
        if roughRT_list[i] != "Missing":
            RTbounds[0].append(roughRT_list[i]-RTtol)
            RTbounds[1].append(roughRT_list[i]+RTtol)
    
    MS1=open(MS1path)
    MS1content=reader(MS1, delimiter=" ")
    MS1filtered=[]
    scan_start = 5
    collect = "no"
    for row in MS1content:
        row = list(row)
        if "S\t" in row[0]: 
            scan_start = int(MS1content.line_num) 
            collect = "no"
        elif "RTime" in row[0]:
            RT_minutes = float(row[0][8:])
            for i in range(0, len(RTbounds[0])):
                if RTbounds[0][i] <= RT_minutes <= RTbounds[1][i]:
                    MS1filtered.append(row)
                    collect = "yes"
                    break       
        elif (MS1content.line_num > (scan_start + 5)) and collect == "yes":
                MS1filtered.append(row)    
    MS1.close()
    
    return MS1filtered

##################################################################################################################################################################

''' 
Generates data for an extracted ion chromatogram (EIC) for a precursor.
Using this data, determines a precise RT by finding the highest point in the chromatogram.
'''

def extraction(pre_mz, pre_mz_tol, MS1list, roughRT, RTtol):

    mz_min = pre_mz - pre_mz/1000000 * pre_mz_tol
    mz_max = pre_mz + pre_mz/1000000 * pre_mz_tol
    EICx, EICy = [] , []
    mzs, intensities = [], []
    RT_minutes = 0
    screen = "no"
    for i in range(0, len(MS1list)):
        if "RTime" in MS1list[i][0]: 
            if i > 0 and screen == "yes":
                if len(mzs) > 1:
                    difference = []
                    for j in range(0, len(mzs)):
                        difference.append(abs(mzs[j] - pre_mz))
                    best_mz_match = difference.index(min(difference))
                    intensity = intensities[best_mz_match]
                    EICx.append(RT_minutes)
                    EICy.append(intensity)
                elif len(mzs) == 1: 
                    intensity = intensities[0]
                    EICx.append(RT_minutes)
                    EICy.append(intensity)  
                else:
                    EICx.append(RT_minutes)
                    EICy.append(0) 
                mzs, intensities = [], []
                screen = "no"
            RT_minutes = float(MS1list[i][0][8:])
            if (roughRT-RTtol) <= RT_minutes <= (roughRT+RTtol):
                screen = "yes"
            elif RT_minutes > (roughRT+RTtol):
                break
        elif screen == "yes": 
            if mz_min <= float(MS1list[i][0]) <= mz_max:
                mzs.append(float(MS1list[i][0]))
                intensities.append(float(MS1list[i][1]))
    peak = max(EICy)
    RT = EICx[EICy.index(peak)]   

    return(RT, EICx, EICy)
    
##################################################################################################################################################################

'''
Using a linear spline approach, generates model to describe relationship between RT in the two runs.
Assuming validation peptide and biological peptide are the same, uses model to predict when validation peptide should elute
    and compares this to actual elution time for validation peptide.
Repeats this for each internal standard (recalculates model with that standard excluded) to determine the accuracy of the model.
'''

def RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT):
    
    RTs = []
    for i in range(0, len(ref_RTs)):
        RTs.append([ref_RTs[i], test_RTs[i]])
        
    sorted_RTs = sorted(RTs, key=lambda entry: entry[0])

    RT_pred_deltas, query_RT_pred_delta = [], 0
    
    for i in range(1, len(sorted_RTs)-1):
        x1, x2 = sorted_RTs[i-1][0], sorted_RTs[i+1][0]
        y1, y2 = sorted_RTs[i-1][1], sorted_RTs[i+1][1] 
        ref, test = sorted_RTs[i][0], sorted_RTs[i][1]
        slope = (y2 - y1)/(x2 - x1)
        y_intercept = y1 - slope*x1
        prediction = slope*ref + y_intercept
        delta = prediction - test
        RT_pred_deltas.append(delta) 
    
    for i in range(0, len(sorted_RTs)):
        if sorted_RTs[i][0] > query_ref_RT:
            hi_index = i
            lo_index = i-1
            break
            
    x1, x2 = sorted_RTs[lo_index][0], sorted_RTs[hi_index][0]
    y1, y2 = sorted_RTs[lo_index][1], sorted_RTs[hi_index][1] 
    ref, test = query_ref_RT, query_test_RT
    slope = (y2 - y1)/(x2 - x1)
    y_intercept = y1 - slope*x1
    prediction = slope*ref + y_intercept
    query_RT_pred_delta = prediction - test
    
    return(sorted_RTs, RT_pred_deltas, query_RT_pred_delta)

##################################################################################################################################################################

'''
ALTERNATIVE RETENTION TIME MODELING APPROACHES...

#cubic spline version

def RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT):
    
    import scipy.interpolate as si
    import copy
    from numpy import ndarray
 
    unique_ref_RTs, unique_test_RTs = [], []
    for i in range(0, len(ref_RTs)):
        if ref_RTs[i] not in unique_ref_RTs:
            unique_ref_RTs.append(ref_RTs[i])
            unique_test_RTs.append(test_RTs[i])
    
    RTs = []
    for i in range(0, len(unique_ref_RTs)):
        RTs.append([unique_ref_RTs[i], unique_test_RTs[i]])
        
    sorted_RTs = sorted(RTs, key=lambda entry: entry[0])

    x, y = [], []
    for i in range(0, len(sorted_RTs)):
        x.append(sorted_RTs[i][0])
        y.append(sorted_RTs[i][1])
        
    query_cs = si.CubicSpline(x, y)
    query_RT_pred = ndarray.tolist(query_cs(query_ref_RT))
    query_RT_pred_delta = query_RT_pred - query_test_RT
    
    RT_pred_deltas = []
    for i in range(1, len(x) - 1):
        ref = copy.deepcopy(x)
        test = copy.deepcopy(y)
        trial_ref = ref.pop(i)
        trial_test = test.pop(i)
        cs = si.CubicSpline(ref, test)
        trial_pred = ndarray.tolist(cs(trial_ref))
        trial_delta = trial_pred - trial_test
        RT_pred_deltas.append(trial_delta)
    
    return(sorted_RTs, RT_pred_deltas, query_RT_pred_delta)    
    
#################################################################################################################################################################

#multi-point version

def RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT):
    
    RTs = []
    for i in range(0, len(ref_RTs)):
        RTs.append([ref_RTs[i], test_RTs[i]])
        
    sorted_RTs = sorted(RTs, key=lambda entry: entry[0])

    RT_pred_deltas, query_RT_pred_delta = [], 0
    
    for i in range(2, len(sorted_RTs)-2):
        x1, x2, x3, x4 = sorted_RTs[i-2][0], sorted_RTs[i-1][0], sorted_RTs[i+1][0], sorted_RTs[i+2][0]
        y1, y2, y3, y4 = sorted_RTs[i-2][1], sorted_RTs[i-1][1], sorted_RTs[i+1][1], sorted_RTs[i+2][1] 
        x_values = [x1, x2, x3, x4]
        y_values = [y1, y2, y3, y4]
        ref, test = sorted_RTs[i][0], sorted_RTs[i][1]
        
        xsquared, xy = [], []
        for i in range(0, len(x_values)):
            xsquared.append(x_values[i]*x_values[i])
            xy.append(x_values[i]*y_values[i])
        sumx = sum(x_values)
        sumy = sum(y_values)
        length = len(x_values)
        y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
        slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
             
        prediction = slope*ref + y_intercept
        delta = prediction - test
        RT_pred_deltas.append(delta) 

    for i in range(0, len(sorted_RTs)):
        if sorted_RTs[i][0] > query_ref_RT:
            hi_index = i
            lo_index = i-1
            break
            
    x1, x2, x3, x4 = sorted_RTs[lo_index - 1][0], sorted_RTs[lo_index][0], sorted_RTs[hi_index][0], sorted_RTs[hi_index + 1][0]
    y1, y2, y3, y4 = sorted_RTs[lo_index - 1][1], sorted_RTs[lo_index][1], sorted_RTs[hi_index][1], sorted_RTs[hi_index + 1][1] 
    x_values = [x1, x2, x3, x4]
    y_values = [y1, y2, y3, y4]
    ref, test = query_ref_RT, query_test_RT
    
    xsquared, xy = [], []
    for i in range(0, len(x_values)):
        xsquared.append(x_values[i]*x_values[i])
        xy.append(x_values[i]*y_values[i])
    sumx = sum(x_values)
    sumy = sum(y_values)
    length = len(x_values)
    y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
    slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
    
    prediction = slope*ref + y_intercept
    query_RT_pred_delta = prediction - test
    
    sorted_RTs_truncated = []
    for i in range(1, len(sorted_RTs)-1):
        sorted_RTs_truncated.append(sorted_RTs[i]) #did this because in main it only removes first and last but for multipoint needs to remove two on each end
    
    return(sorted_RTs_truncated, RT_pred_deltas, query_RT_pred_delta)

#################################################################################################################################################################

#basic linear model
    
import copy

def RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT):

    RTs = []
    for i in range(0, len(ref_RTs)):
        RTs.append([ref_RTs[i], test_RTs[i]])
        
    sorted_RTs = sorted(RTs, key=lambda entry: entry[0])
    
    sorted_ref_RTs, sorted_test_RTs = [], []
    for i in range(0, len(sorted_RTs)):
        sorted_ref_RTs.append(sorted_RTs[i][0])
        sorted_test_RTs.append(sorted_RTs[i][1])
    
    RT_pred_deltas, query_RT_pred_delta = [], 0
    
    for i in range(1, len(sorted_ref_RTs)-1):
        x_values = copy.deepcopy(sorted_ref_RTs)
        y_values = copy.deepcopy(sorted_test_RTs)

        ref = x_values.pop(i)
        test = y_values.pop(i)

        xsquared, xy = [], []
        for i in range(0, len(x_values)):
            xsquared.append(x_values[i]*x_values[i])
            xy.append(x_values[i]*y_values[i])
        sumx = sum(x_values)
        sumy = sum(y_values)
        length = len(x_values)
        y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
        slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
             
        prediction = slope*ref + y_intercept
        delta = prediction - test
        RT_pred_deltas.append(delta) 
        
    x_values = sorted_ref_RTs
    y_values = sorted_test_RTs

    xsquared, xy = [], []
    for i in range(0, len(x_values)):
        xsquared.append(x_values[i]*x_values[i])
        xy.append(x_values[i]*y_values[i])
    sumx = sum(x_values)
    sumy = sum(y_values)
    length = len(x_values)
    y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
    slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
         
    prediction = slope*query_ref_RT + y_intercept
    query_RT_pred_delta = prediction - query_test_RT
    
    return(sorted_RTs, RT_pred_deltas, query_RT_pred_delta)

#################################################################################################################################################################

#log2 linear model
    
import copy

def RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT):

    import math
    
    RTs = []
    for i in range(0, len(ref_RTs)):
        RTs.append([ref_RTs[i], test_RTs[i]])
        
    sorted_RTs = sorted(RTs, key=lambda entry: entry[0])
    
    sorted_ref_RTs, sorted_test_RTs = [], []
    for i in range(0, len(sorted_RTs)):
        sorted_ref_RTs.append(math.log(sorted_RTs[i][0], 2))
        sorted_test_RTs.append(math.log(sorted_RTs[i][1], 2))
    
    RT_pred_deltas, query_RT_pred_delta = [], 0
    
    for i in range(1, len(sorted_ref_RTs)-1):
        x_values = copy.deepcopy(sorted_ref_RTs)
        y_values = copy.deepcopy(sorted_test_RTs)

        ref = x_values.pop(i)
        test = y_values.pop(i)

        xsquared, xy = [], []
        for i in range(0, len(x_values)):
            xsquared.append(x_values[i]*x_values[i])
            xy.append(x_values[i]*y_values[i])
        sumx = sum(x_values)
        sumy = sum(y_values)
        length = len(x_values)
        y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
        slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
             
        prediction = slope*ref + y_intercept
        delta = prediction - test
        RT_pred_deltas.append(delta) 
        
    x_values = sorted_ref_RTs
    y_values = sorted_test_RTs

    xsquared, xy = [], []
    for i in range(0, len(x_values)):
        xsquared.append(x_values[i]*x_values[i])
        xy.append(x_values[i]*y_values[i])
    sumx = sum(x_values)
    sumy = sum(y_values)
    length = len(x_values)
    y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
    slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
         
    prediction = slope*math.log(query_ref_RT, 2) + y_intercept
    query_RT_pred_delta = prediction - math.log(query_test_RT, 2)
    
    return(sorted_RTs, RT_pred_deltas, query_RT_pred_delta)

#################################################################################################################################################################

#log2 linear model with back-transform
    
import copy

def RT_prediction(ref_RTs, test_RTs, query_ref_RT, query_test_RT):
    
    import math
        

    RTs = []
    for i in range(0, len(ref_RTs)):
        RTs.append([ref_RTs[i], test_RTs[i]])
        
    sorted_RTs = sorted(RTs, key=lambda entry: entry[0])
    
    sorted_ref_RTs, sorted_test_RTs = [], []
    for i in range(0, len(sorted_RTs)):
        sorted_ref_RTs.append(math.log(sorted_RTs[i][0], 2))
        sorted_test_RTs.append(math.log(sorted_RTs[i][1], 2))
    
    RT_pred_deltas, query_RT_pred_delta = [], 0
    
    for i in range(1, len(sorted_ref_RTs)-1):
        x_values = copy.deepcopy(sorted_ref_RTs)
        y_values = copy.deepcopy(sorted_test_RTs)

        ref = x_values.pop(i)
        test = y_values.pop(i)

        xsquared, xy = [], []
        for i in range(0, len(x_values)):
            xsquared.append(x_values[i]*x_values[i])
            xy.append(x_values[i]*y_values[i])
        sumx = sum(x_values)
        sumy = sum(y_values)
        length = len(x_values)
        y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
        slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
             
        prediction = slope*ref + y_intercept
        delta = 2**prediction - 2**test
        RT_pred_deltas.append(delta) 
        
    x_values = sorted_ref_RTs
    y_values = sorted_test_RTs

    xsquared, xy = [], []
    for i in range(0, len(x_values)):
        xsquared.append(x_values[i]*x_values[i])
        xy.append(x_values[i]*y_values[i])
    sumx = sum(x_values)
    sumy = sum(y_values)
    length = len(x_values)
    y_intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
    slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
         
    prediction = slope*math.log(query_ref_RT, 2) + y_intercept
    query_RT_pred_delta = 2**prediction - query_test_RT
    
    return(sorted_RTs, RT_pred_deltas, query_RT_pred_delta)

#################################################################################################################################################################

'''