from math import atanh, sqrt, tanh
from statistics import stdev
from statistics import mean as statmean
from scipy.stats import norm, normaltest
from scipy.stats import t as stat_t
import matplotlib.pyplot as plt

##################################################################################################################################################################

def PCC_percentile_rank(r_list, percentile_thresh, query_r, out_dir, file_name):
    
#r_list: list of Pearson r values for internal standards
#conf_level: user-defined level for prediction interval (e.g., 0.9, 0.95)    
#percentile_thresh: user-defined minimum percentile rank; the Pearson r of the query must exceed this rank to pass 
#query_r: Pearson r for the peptide of interest
    
    #Fisher's z transformation of Pearson r data
    
    z_list = []
    for i in range(0, len(r_list)):
            #The traditional equation is z = 0.5 * (math.log(1+data[i]) - math.log(1-data[i])), but it can also be determined by calculating the inverse
            #hyperbolic tangent (which makes for cleaner code). I confirmed that both approaches give the same answer.
        z = atanh(r_list[i])
        z_list.append(z)
        
    #Calculate summary statistics
    
    mean_z = statmean(z_list)
    std_z = stdev(z_list)
    n_z = len(z_list)

    #Test normality
    
#    alpha = 0.05
    normality = normaltest(z_list)
    normality_p = float(normality[1])
#    if normality_p >= alpha:
#        normality_check = "passed normality test"
#    else:
#        normality_check = "failed normality test"
        
    plt.figure(figsize=(5,5))
    plt.hist(z_list, bins=10, density=True)
    z_range = max(z_list) - min(z_list)
    step = z_range/100
    x = [min(z_list)]
    for i in range(1,101):
        x.append(x[i-1]+step)
    plt.plot(x, norm(mean_z,std_z).pdf(x))
    plt.xlabel("Fisher z", fontsize="x-large")
    plt.ylabel("probability density", fontsize="x-large")
    plt.xticks(fontsize="large")
    plt.yticks(fontsize="large")
    plt.title("normality of ISP PCCs (p = "+str(round(normality_p,3))+")", fontsize="xx-large", pad = 15)
    plt.savefig(out_dir + "\\Figures\\normality\\" + file_name + "_PCC_normality", bbox_inches = "tight") 
    plt.close(fig='all')
    
    #Calculate t-statistic for threshold
    
    t_thresh = stat_t.ppf(percentile_thresh, n_z - 1) #One-tailed percentile approach; percentile_thresh = q
    
    z_thresh = mean_z + t_thresh*std_z*sqrt(1 + 1/n_z)
    
    #Back transform z_thresh to obtain minimum Pearson r
    
    r_thresh = tanh(z_thresh)
        
    #Calculate percentile for query_r
    
    query_z = atanh(query_r)
    query_t = (query_z - mean_z)/(std_z*sqrt(1 + 1/n_z))
    query_percentile = stat_t.cdf(query_t, n_z - 1)*100

    return(r_thresh, query_percentile)
    #r_lo and r_hi define the boundaries for the prediction interval

##################################################################################################################################################################

def RT_percentile_rank(RT_pred_deltas, percentile_thresh, query_RT_pred_delta, out_dir, file_name):
    
    #Calculate summary statistics
    
    mean = statmean(RT_pred_deltas)
    std = stdev(RT_pred_deltas)
    n = len(RT_pred_deltas)

    #Test normality
    
#    alpha = 0.05
    normality = normaltest(RT_pred_deltas)
    normality_p = float(normality[1])
#    if normality_p >= alpha:
#        normality_check = "passed normality test"
#    else:
#        normality_check = "failed normality test"
        
    plt.figure(figsize=(5,5))
    plt.hist(RT_pred_deltas, bins=10, density=True)
    delta_range = max(RT_pred_deltas) - min(RT_pred_deltas)
    step = delta_range/100
    x = [min(RT_pred_deltas)]
    for i in range(1,101):
        x.append(x[i-1]+step)
    plt.plot(x, norm(mean,std).pdf(x))
    plt.xlabel("expected RT - syn RT (minutes)", fontsize="x-large")
    plt.ylabel("probability density", fontsize="x-large")
    plt.xticks(fontsize="large")
    plt.yticks(fontsize="large")
    plt.title("normality of ISP delta RTs (p = "+str(round(normality_p,3))+")", fontsize="xx-large", pad = 15)
    plt.savefig(out_dir + "\\Figures\\normality\\" + file_name + "_deltaRT_normality", bbox_inches = "tight") 
    plt.close(fig='all')
    
    #Calculate t-statistic and prediction interval boundaries
    
    t = stat_t.ppf(percentile_thresh/2, n - 1) #This is for two-tailed calculation. PPF: percent point function (quantile function).
    
    delta_lo = mean + t*std*sqrt(1 + 1/n) 
    delta_hi = mean - t*std*sqrt(1 + 1/n) 

    #Calculate percentile rank for query_RT_pred_delta
    
    query_t = (query_RT_pred_delta - mean)/(std*sqrt(1 + 1/n))
    query_percentile = stat_t.cdf(-abs(query_t), n - 1)*100*2 #Multiply by two to provide two-tailed value

    return(delta_lo, delta_hi, query_percentile)

##################################################################################################################################################################