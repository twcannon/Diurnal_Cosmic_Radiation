import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import os,csv,sys
import math
from scipy.stats import spearmanr


### READ GEIGER DATA ###
os.chdir('/Users/Thomas/Desktop/Python')
reader=csv.reader(open("cosmic_ray_data.csv","rU"),delimiter=',')
x=list(reader)
result=np.array(x)
result_matrix = result[ range(1 , len(result[:,0]) )  ,  :].astype(float)

### ARRAY ORIENTATION ###
setup = result_matrix[:,0]                       
date = result_matrix[:,1]                        
time = result_matrix[:,2]                     
counts = result_matrix[:,3]                     
duration = result_matrix[:,4]                        
pressure = result_matrix[:,5]                     
temp = result_matrix[:,6]                     
humidity = result_matrix[:,7]   

#loc_20 = np.where(duration == 20)
#for i in range(0, len(loc_20[0])):
#    counts[i] = counts[i]/4
#    

                     
                                          
 ### READ SUNSPOT DATA ###
os.chdir('/Users/Thomas/Desktop/Python')
reader=csv.reader(open("sunspot_data.csv","rU"),delimiter=',')
x=list(reader)
sunspot=np.array(x)
sunspot_matrix = sunspot[ range(1 , len(sunspot[:,0]) )  ,  :].astype(float)                     

### ARRAY ORIENTATION ###                    
date_sunspot = sunspot_matrix[:,0]    
sunspots = sunspot_matrix[:,1]                       



### DETECTOR ONE ###
loc1 = np.where(setup == 1)
for i in range(len(loc1[0])):
    setup1 = setup[loc1[0:i]]
    date1 = date[loc1[0:i]]
    time1 = time[loc1[0:i]]
    counts1 = counts[loc1[0:i]]
    duration1 = duration[loc1[0:i]]
    pressure1 = pressure[loc1[0:i]]
    temp1 = temp[loc1[0:i]]
    humidity1 = humidity[loc1[0:i]]
beq1 = counts1/(duration1*60.)    
    
    
### DETECTOR TWO ###   
loc2 = np.where(setup == 2) 
for i in range(len(loc2[0])):
    setup2 = setup[loc2[0:i]]
    date2 = date[loc2[0:i]]
    time2 = time[loc2[0:i]]
    counts2 = counts[loc2[0:i]]
    duration2 = duration[loc2[0:i]]
    pressure2 = pressure[loc2[0:i]]
    temp2 = temp[loc2[0:i]]
    humidity2 = humidity[loc2[0:i]]    
beq2 = counts2/(duration2*60.)


### DETECTOR THREE ###   
loc3 = np.where(setup == 3) 
for i in range(len(loc3[0])):
    setup3 = setup[loc3[0:i]]
    date3 = date[loc3[0:i]]
    time3 = time[loc3[0:i]]
    counts3 = counts[loc3[0:i]]
    duration3 = duration[loc3[0:i]]
    pressure3 = pressure[loc3[0:i]]
    temp3 = temp[loc3[0:i]]
    humidity3 = humidity[loc3[0:i]]    
beq3 = counts3/(duration3*60.)


### DETECTOR FOUR ###   
loc4 = np.where(setup == 4) 
for i in range(len(loc4[0])):
    setup4 = setup[loc4[0:i]]
    date4 = date[loc4[0:i]]
    time4 = time[loc4[0:i]]
    counts4 = counts[loc4[0:i]]
    duration4 = duration[loc4[0:i]]
    pressure4 = pressure[loc4[0:i]]
    temp4 = temp[loc4[0:i]]
    humidity4 = humidity[loc4[0:i]]    
beq4 = counts4/(duration4*60.)


beq = np.concatenate((beq1,beq2,beq3,beq4))
#time = np.concatenate((time1,time2,time3,time4))
pressure = np.concatenate((pressure1,pressure2,pressure3,pressure4))
temp = np.concatenate((temp1,temp2,temp3,temp4))
humidity = np.concatenate((humidity1,humidity2,humidity3,humidity4))


### NORMALIZATION ###

beq_av = np.mean(beq)
beq_lab_av = np.mean(np.concatenate((beq1,beq2)))
beq_home_av = np.mean(np.concatenate((beq3,beq4)))
beq1_av = np.mean(beq1) 
beq2_av = np.mean(beq2) 
beq3_av = np.mean(beq3) 
beq4_av = np.mean(beq4)
beq1_n = (beq1/beq1_av)*beq_lab_av 
beq2_n = (beq2/beq2_av)*beq_lab_av 
beq_lab_norm = np.concatenate((beq1_n,beq2_n))
beq3_n = (beq3/beq3_av)*beq_home_av 
beq4_n = (beq4/beq4_av)*beq_home_av 
beq_lab_norm = np.concatenate((beq1_n,beq2_n))
beq_home_norm = np.concatenate((beq3_n,beq4_n))

beq_home_norm_av = np.mean(beq_home_norm)
beq_lab_norm_av = np.mean(beq_lab_norm)
beq1_norm = (beq1_n/beq_lab_norm_av)*beq_av
beq2_norm = (beq2_n/beq_lab_norm_av)*beq_av
beq3_norm = (beq3_n/beq_home_norm_av)*beq_av
beq4_norm = (beq4_n/beq_home_norm_av)*beq_av

#beq1_norm = (beq1_n/beq1_av)*beq_av
#beq2_norm = (beq2_n/beq2_av)*beq_av
#beq3_norm = (beq3_n/beq3_av)*beq_av
#beq4_norm = (beq4_n/beq4_av)*beq_av

beq_norm = np.concatenate((beq1_norm,beq2_norm,beq3_norm,beq4_norm))


loc_day = np.where(np.logical_and(time <=  19.75, time >=  5))
loc_night = np.where(np.logical_or(time <=  5,time >  19.75))
beq_day = beq_norm[loc_day]  
setup_day = setup[loc_day]
date_day = date[loc_day]

loc1 = np.where(setup_day == 1) 
date_1 = date_day[loc1]
loc2 = np.where(setup_day == 2) 
date_2 = date_day[loc2]
loc3 = np.where(setup_day == 3) 
date_3 = date_day[loc3]
loc4 = np.where(setup_day == 4) 
date_4 = date_day[loc4]

dates = np.concatenate((date_1,date_2,date_3,date_4))

means = []
for i in range(0,len(date_sunspot)):
    loc = np.where(dates ==  date_sunspot[i])
    mean_beq_date = np.mean(beq_norm[loc])
    means.append(mean_beq_date)
    
xerr_arr = np.sqrt(sunspots)
yerr_arr = []
for i in range(0,len(means)):
    y_err = np.sqrt(means[i]*5.*60.)/(5.*60.)
    yerr_arr.append(y_err)


rho_p, P_p = spearmanr(means,sunspots)
print 'rho:',rho_p
print 'probability:',round((1-P_p)*100,1),'%'


plt.plot(sunspots,means,'bs')
#plt.errorbar(sunspots, means, yerr=yerr_arr, xerr=xerr_arr, fmt='bs')#, ecolor=None, elinewidth=None, capsize=3)
plt.xlabel('Number of Sunspots')
plt.ylabel('Becquerel')
plt.ylim(.2,.5)
plt.show()