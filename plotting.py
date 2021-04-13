import matplotlib.pyplot as plt
from copy import deepcopy

##################################################################################################################################################################
    
def swarm_plot(y_list, series_name, y_highlight, highlight_series_name, title, y_axis_label, out_dir, file_name, file_suffix, lo_thresh, hi_thresh, percentile):
    
    y_max = max(max(y_list), y_highlight, lo_thresh, hi_thresh)
    y_min = min(min(y_list), y_highlight, lo_thresh, hi_thresh)
    y_range = y_max - y_min

    x_list=[]
    for i in range(0,len(y_list)):
        y=y_list[i]
        similarity=0
        for j in range (0,i):
            if abs(y_list[j] - y) < 0.03*y_range:
                similarity = similarity + 1
        
        if similarity > 0:
            if similarity%2 == 0:
                x_list.append(1+similarity/2*.05) 
            elif similarity%2 == 1:
                x_list.append(1-(similarity+1)/2*.05)
        else:
            x_list.append(1)
    
    plt.figure(figsize=(5,5))
    plt.fill_between([0.5*min(x_list), 1.5*max(x_list)], lo_thresh, hi_thresh, facecolor="green", alpha=.1)
    plt.scatter(x_list,y_list,label= series_name, color = "black")
    plt.scatter(1, y_highlight, label= highlight_series_name, color = "red", marker = "s", s=60)
    plt.xlim(min(x_list)-.4,max(x_list)+.4)
    plt.ylim(y_min - 0.3*y_range, y_max + 0.1*y_range)
    plt.text(max(x_list)+0.375,y_max,"%tile="+str(percentile)+"%",fontsize="x-large", horizontalalignment="right")     
    plt.xticks([])
    plt.yticks(fontsize="large")
    plt.legend(loc="lower left", fontsize="large")   
    plt.title(title, fontsize="xx-large", pad = 15)
    plt.ylabel(y_axis_label, fontsize="x-large")
    plt.savefig(out_dir + file_name + file_suffix, bbox_inches = "tight") 
    plt.close(fig='all')
    
##################################################################################################################################################################

def regression_plot(x_values, y_values, series_name, query_x, query_y, query_name, title, x_label, y_label, regression, scientific, text, out_dir, filename, vline1, vline2):

    plt.figure(figsize=(5,5))
    plt.scatter(x_values,y_values, label = series_name, color="black")
    if query_x != "NA":
        y_max = max(max(y_values), query_y)
        y_min = min(min(y_values), query_y)
    else:
        y_max = max(y_values)
        y_min = min(y_values)
    y_range = y_max - y_min
    y_bottom = y_min - 0.3*y_range
    y_top = y_max + 0.1*y_range
    if vline1 != "N" and vline2 != "N":
        plt.vlines([vline1, vline2], [y_bottom, y_bottom], [y_top, y_top], linestyles="dotted", lw=2, color = "black")
    if query_x != "NA":
        plt.scatter(query_x, query_y, label = query_name, color = "red", marker = "s", s = 60)
        plt.ylim(y_bottom, y_top)
        plt.legend(loc="lower right", fontsize="large")
    plt.title(title, fontsize="xx-large", pad = 15)
    plt.xlabel(x_label, fontsize="x-large")
    plt.ylabel(y_label, fontsize="x-large")
    plt.xticks(fontsize="large")
    plt.yticks(fontsize="large")
    if scientific == "Y":
        plt.ticklabel_format(style="sci", scilimits=(0,0), axis="both")
    #regression line
    if regression == "Y":
        xsquared, xy = [], []
        for i in range(0, len(x_values)):
            xsquared.append(x_values[i]*x_values[i])
            xy.append(x_values[i]*y_values[i])
        sumx = sum(x_values)
        sumy = sum(y_values)
        length = len(x_values)
        intercept = (sumy*sum(xsquared)-sumx*sum(xy))/(length*sum(xsquared)-sumx*sumx)
        slope =  (length*sum(xy)-sumx*sumy)/(length*sum(xsquared) - sumx*sumx)
        fitx = [min(x_values), max(x_values)]
        fity = [slope*min(x_values) + intercept, slope*max(x_values) + intercept]
        plt.plot(fitx, fity, marker = "", color = "black", linestyle = "--")
    if text != "N":
        plt.text(max(x_values), min(y_values),text,fontsize="x-large", horizontalalignment="right")   
    plt.savefig(out_dir + filename + ".png", bbox_inches = "tight")
    plt.close(fig='all')
    
##################################################################################################################################################################

def mirror_plot(hyphen, sequence, N_term_shift, C_term_shift, abund_thresh, biomin, synmin, bio_scan, bio_L_ions, bio_R_ions, syn_scan, syn_L_ions, syn_R_ions, ion_type, out_dir, filename):
    
    biox,bioy=[],[]
    bioLx,bioLy=[],[]
    synx,syny=[],[]
    synLx,synLy=[],[]
    xlimits=[]
    if hyphen != -1:
        bioRx,bioRy=[],[]
        synRx,synRy=[],[]
    
    for i in range(5, len(bio_scan)-1):
        biox.append(float(bio_scan[i][0]))
        bioy.append(float(bio_scan[i][1]))
    bioy_max=max(bioy)
    bioy_norm=deepcopy(bioy)
    for i in range(len(bioy_norm)):
        bioy_norm[i]=bioy_norm[i]/bioy_max
    for i in range(len(biox)):
        if bioy[i] >= biomin and bioy[i] >= abund_thresh:
            biox_max = biox[i]
    xlimits.append(biox_max)        
    xlimits.append(min(biox))
        
    for i in range(5, len(syn_scan)-1):
        synx.append(float(syn_scan[i][0]))
        syny.append(float(syn_scan[i][1])) 
    syny_max=max(syny)
    syny_norm=deepcopy(syny)
    for i in range(len(syny_norm)):
        syny_norm[i]=-syny_norm[i]/syny_max
    for i in range(len(synx)):
        if syny[i] >= synmin and syny[i] >= abund_thresh:
            synx_max = synx[i]
    xlimits.append(synx_max)       
    xlimits.append(min(synx))
    
    LRmap=[]
    
    for i in range(len(bio_L_ions)):
        b, y = "", ""
        if bio_L_ions[i][7] != "":
            bioLx.append(float(bio_L_ions[i][7]))
            bioLy.append(float(bio_L_ions[i][12])/bioy_max)
            b = "found"
        if bio_L_ions[i][8] != "":
            bioLx.append(float(bio_L_ions[i][8]))
            bioLy.append(float(bio_L_ions[i][13])/bioy_max)
            b = "found"
        if bio_L_ions[i][9] != "":
            bioLx.append(float(bio_L_ions[i][9]))
            bioLy.append(float(bio_L_ions[i][14])/bioy_max)
            y = "found"
        if bio_L_ions[i][10] != "":
            bioLx.append(float(bio_L_ions[i][10]))
            bioLy.append(float(bio_L_ions[i][15])/bioy_max)
            y = "found"
        if b == "found" and y == "found":
            LRmap.append("both")
        elif b == "found":
            LRmap.append("b")
        elif y == "found":
            LRmap.append("y")
        else:
            LRmap.append("none")
    if hyphen != -1:
        for i in range(len(bio_R_ions)):
            b, y = "", ""
            if bio_R_ions[i][7] != "":
                bioRx.append(float(bio_R_ions[i][7]))
                bioRy.append(float(bio_R_ions[i][12])/bioy_max)
                b = "found"
            if bio_R_ions[i][8] != "":
                bioRx.append(float(bio_R_ions[i][8]))
                bioRy.append(float(bio_R_ions[i][13])/bioy_max)
                b = "found"
            if bio_R_ions[i][9] != "":
                bioRx.append(float(bio_R_ions[i][9]))
                bioRy.append(float(bio_R_ions[i][14])/bioy_max)
                y = "found"
            if bio_R_ions[i][10] != "":
                bioRx.append(float(bio_R_ions[i][10]))
                bioRy.append(float(bio_R_ions[i][15])/bioy_max)
                y = "found"
            if b == "found" and y == "found":
                LRmap.append("both")
            elif b == "found":
                LRmap.append("b")
            elif y == "found":
                LRmap.append("y")
            else:
                LRmap.append("none")
                
    sequence_map=sequence[0]
    for i in range(1, len(sequence)):
        if LRmap[i-1] == "both":
            sequence_map=sequence_map+"|"
        elif LRmap[i-1] == "b":
            sequence_map=sequence_map+"\\"
        elif LRmap[i-1] == "y":
            sequence_map=sequence_map+"/"
        elif LRmap[i-1] == "none":
            sequence_map=sequence_map+" "
        sequence_map=sequence_map+sequence[i]
    if N_term_shift != 0:
        sequence_map = "("+str(round(N_term_shift,2))+")"+sequence_map 
    if C_term_shift != 0:
        sequence_map = sequence_map+"("+str(round(C_term_shift,2))+")"    
            
    for i in range(len(syn_L_ions)):
        if syn_L_ions[i][7] != "":
            synLx.append(float(syn_L_ions[i][7]))
            synLy.append(-float(syn_L_ions[i][12])/syny_max)
        if syn_L_ions[i][8] != "":
            synLx.append(float(syn_L_ions[i][8]))
            synLy.append(-float(syn_L_ions[i][13])/syny_max)
        if syn_L_ions[i][9] != "":
            synLx.append(float(syn_L_ions[i][9]))
            synLy.append(-float(syn_L_ions[i][14])/syny_max)
        if syn_L_ions[i][10] != "":
            synLx.append(float(syn_L_ions[i][10]))
            synLy.append(-float(syn_L_ions[i][15])/syny_max)
    if hyphen != -1:
        for i in range(len(syn_R_ions)):
            if syn_R_ions[i][7] != "":
                synRx.append(float(syn_R_ions[i][7]))
                synRy.append(-float(syn_R_ions[i][12])/syny_max)
            if syn_R_ions[i][8] != "":
                synRx.append(float(syn_R_ions[i][8]))
                synRy.append(-float(syn_R_ions[i][13])/syny_max)
            if syn_R_ions[i][9] != "":
                synRx.append(float(syn_R_ions[i][9]))
                synRy.append(-float(syn_R_ions[i][14])/syny_max)
            if syn_R_ions[i][10] != "":
                synRx.append(float(syn_R_ions[i][10]))
                synRy.append(-float(syn_R_ions[i][15])/syny_max)
    
    plt.figure(figsize=(10,5))
    plt.yticks([])
    plt.xticks(fontsize = "large")
    #plt.xticks(fontsize = 20, rotation = 90, fontweight = "bold") #instead of bold can do a numeric value in range 0-1000, 'ultralight', 'light', 'normal', 'regular', 'book', 'medium', 'roman', 'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'
    plt.title(sequence_map, fontsize=21, fontweight = "semibold", fontname="Courier New", pad = 12) #Consolas; I like Comic Sans MS
    plt.ylabel("normalized spectral intensity", fontsize="x-large")
    plt.xlabel("m/z", fontsize="x-large")
    plt.stem(biox,bioy_norm,"k", markerfmt=" ", basefmt="k-", use_line_collection=True)
    plt.stem(synx,syny_norm,"k", markerfmt=" ", basefmt="k-", use_line_collection=True)
    stemplot = plt.stem(bioLx,bioLy,"r", markerfmt="ro", basefmt="k-", use_line_collection=True)
    marker = stemplot[0]
    plt.setp(marker, markersize = 5)
    if hyphen == -1:
        if ion_type == "b/y":
            stemplot = plt.stem(synLx,synLy,"r", markerfmt="ro", basefmt="k-", label="b/y ions", use_line_collection=True)
        elif ion_type == "c/z":
            stemplot = plt.stem(synLx,synLy,"r", markerfmt="ro", basefmt="k-", label="c/z ions", use_line_collection=True)            
    else:
        stemplot = plt.stem(synLx,synLy,"r", markerfmt="ro", basefmt="k-", label="L ions", use_line_collection=True)        
    marker = stemplot[0]
    plt.setp(marker, markersize = 5)
    if hyphen != -1:
        stemplot = plt.stem(bioRx,bioRy,"c", markerfmt="co", basefmt="k-", use_line_collection=True)
        marker = stemplot[0]
        plt.setp(marker, markersize = 5)
        stemplot = plt.stem(synRx,synRy,"c", markerfmt="co", basefmt="k-", label="R ions", use_line_collection=True)
        marker = stemplot[0]
        plt.setp(marker, markersize = 5)        
    plt.legend(loc="lower right", fontsize="x-large")
    bio_hline = max(biomin/bioy_max, abund_thresh/bioy_max)
    syn_hline = max(synmin/syny_max, abund_thresh/syny_max)
    plt.hlines(y=bio_hline, xmin=min(xlimits), xmax=max(xlimits), linestyles="dotted", lw=3, color = "black")
    plt.hlines(y=-syn_hline, xmin=min(xlimits), xmax=max(xlimits), linestyles="dotted", lw=3, color = "black")
    xrange = max(xlimits) - min(xlimits)
    if hyphen == -1:
        plt.text(max(xlimits) + 0.025*xrange,-.65,f"synthetic\n(max={syny_max:.1E})", fontsize="x-large", horizontalalignment="right")
    else:
        plt.text(max(xlimits) + 0.025*xrange,-.55,f"synthetic\n(max={syny_max:.1E})", fontsize="x-large", horizontalalignment="right")        
    plt.text(max(xlimits) + 0.025*xrange,.8,f"biological\n(max={bioy_max:.1E})", fontsize="x-large", horizontalalignment="right")
    plt.xlim(min(xlimits) - 0.05*xrange, max(xlimits) + 0.05*xrange)
    plt.savefig(out_dir + filename, bbox_inches = "tight")
    plt.close(fig='all')
    
    return (sequence_map, biox, bioy_max, bioy_norm, synx, syny_max, syny_norm)

##################################################################################################################################################################

def EIC(series1_x, series1_y, series1_RT, series1_label, series2_x, series2_y, series2_RT, series2_label, compoundID, out_dir, filename):
    
    plt.figure(figsize=(8,5))  
    fig,ax = plt.subplots()
    
    ax.plot(series1_x, series1_y, color="blue")
    ax.set_xlabel("retention time (minutes)", fontsize="x-large")
    ax.set_ylabel("raw intensity in "+series1_label,color="blue",fontsize="x-large")
    plt.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    plt.xticks(fontsize="large")
    plt.yticks(fontsize="large")
    plt.vlines(x = series1_RT, ymin = 0, ymax = max(series1_y), linestyles = "dashed", color = "black")
    plt.text(series1_RT, 0.03*max(series1_y), round(series1_RT, 1), color = "blue", fontsize = "x-large", horizontalalignment = "center", bbox = dict(facecolor = "white", edgecolor = "None"))
    
    ax2=ax.twinx()
    ax2.plot(series2_x, series2_y,color="red")
    ax2.set_ylabel("raw intensity in "+series2_label,color="red",fontsize="x-large")
    plt.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    plt.xticks(fontsize="large")
    plt.yticks(fontsize="large")
    plt.vlines(x = series2_RT, ymin = 0, ymax = max(series2_y), linestyles = "dashed", color = "black")
    plt.text(series2_RT, 0.03*max(series2_y), round(series2_RT, 1), color = "red", fontsize = "x-large", horizontalalignment = "center", bbox = dict(facecolor = "white", edgecolor = "None"))
    
    plt.title("extracted ion chromatograms: " + compoundID, fontsize="xx-large", pad = 15)  

    plt.savefig(out_dir + filename, bbox_inches = "tight") 
    plt.close(fig='all')
    
##################################################################################################################################################################
