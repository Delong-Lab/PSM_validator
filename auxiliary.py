##################################################################################################################################################################

def time_format(now):

    year = str(now.year)
    year = year[2:]
    month = str(now.month)
    if len(month) == 1:
        month = "0" + month
    day = str(now.day)
    if len(day) == 1:
        day = "0" + day
    hour = str(now.hour)
    if len(hour) == 1:
        hour = "0" + hour
    minute = str(now.minute)
    if len(minute) == 1:
        minute = "0" + minute
    second = str(now.second)
    if len(second) == 1:
        second = "0" + second
    
    date = year + month + day 
    time = hour + minute + second
    
    return (date, time)

##################################################################################################################################################################

def timer(now, start_time):

    processing_time = now-start_time
    minutes = str(int(processing_time.seconds)//60)
    seconds = str(int(processing_time.seconds)%60)
    microseconds = str(processing_time.microseconds)
    processing_time = minutes+":"+seconds+"."+microseconds
    
    return (processing_time)

##################################################################################################################################################################