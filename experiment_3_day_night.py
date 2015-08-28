import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import os,csv,sys
import math


### READ DATA ###
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

loc_24_beq = np.where(date ==  24)
beq_24 = beq_norm[loc_24_beq]    
loc_26_beq = np.where(date ==  26)
beq_26 = beq_norm[loc_26_beq]  
loc_32_beq = np.where(date ==  32)
beq_32 = beq_norm[loc_32_beq] 
    
loc_24_3_beq = np.where(setup[loc_24_beq] ==  3)
beq_24_3 = beq_24[loc_24_3_beq]  
loc_24_4_beq = np.where(setup[loc_24_beq] ==  4)
beq_24_4 = beq_24[loc_24_4_beq]     
loc_26_3_beq = np.where(setup[loc_26_beq] ==  3)
beq_26_3 = beq_26[loc_26_3_beq]  
loc_26_4_beq = np.where(setup[loc_26_beq] ==  4)
beq_26_4 = beq_26[loc_26_4_beq] 
loc_32_3_beq = np.where(setup[loc_32_beq] ==  3)
beq_32_3 = beq_32[loc_32_3_beq]  
loc_32_4_beq = np.where(setup[loc_32_beq] ==  4)
beq_32_4 = beq_32[loc_32_4_beq] 
    
    
loc_24_time = np.where(date ==  24)
time_24 = time[loc_24_time]    
loc_26_time = np.where(date ==  26)
time_26 = time[loc_26_time] 
loc_32_time = np.where(date ==  32)
time_32 = time[loc_32_time] 
    
    
loc_24_3_time = np.where(setup[loc_24_time] ==  3)
time_24_3 = time_24[loc_24_3_time]     
loc_24_4_time = np.where(setup[loc_24_time] ==  4)
time_24_4 = time_24[loc_24_4_time]     
loc_26_3_time = np.where(setup[loc_26_time] ==  3)
time_26_3 = time_26[loc_26_3_time]     
loc_26_4_time = np.where(setup[loc_26_time] ==  4)
time_26_4 = time_26[loc_26_4_time]   
loc_32_3_time = np.where(setup[loc_32_time] ==  3)
time_32_3 = time_32[loc_32_3_time]     
loc_32_4_time = np.where(setup[loc_32_time] ==  4)
time_32_4 = time_32[loc_32_4_time]   
    

med_beq = []
#for i in range(0,len(time_26_4)):
#    med_array = timenp.where()
    
for i in range(0,1):
    med_array = np.array((beq_26_3[i],beq_26_4[i]))
    med_beq.append(np.median(med_array))
for i in range(1,2):
    med_array = np.array((beq_26_3[i],beq_26_4[i],beq_32_3[i-1],beq_32_4[i-1]))
    med_beq.append(np.median(med_array))
for i in range(0,len(time_24_4)):
    med_array = np.array((beq_24_3[i],beq_24_4[i],beq_26_3[i+2],beq_26_4[i+2],beq_32_3[i+1],beq_32_4[i+1]))
    med_beq.append(np.median(med_array))
    
mean_beq = []
for i in range(0,2):
    mean_array = np.array((beq_26_3[i],beq_26_4[i],beq_32_3[i],beq_32_4[i]))
    mean_beq.append(np.mean(mean_array))
for i in range(1,2):
    mean_array = np.array((beq_26_3[i],beq_26_4[i],beq_32_3[i-1],beq_32_4[i-1]))
    mean_beq.append(np.median(mean_array))
for i in range(1,len(time_24_4)):
    mean_array = np.array((beq_24_3[i],beq_24_4[i],beq_26_3[i+2],beq_26_4[i+2],beq_32_3[i+1],beq_32_4[i+1]))
    mean_beq.append(np.mean(mean_array))    


plt.plot(time_24_3,beq_24_3, color = 'c')
plt.plot(time_24_4,beq_24_4, color = 'm')
plt.plot(time_26_3,beq_26_3, color = 'c')
plt.plot(time_26_4,beq_26_4, color = 'm')
plt.plot(time_32_3,beq_32_3, color = 'c')
plt.plot(time_32_4,beq_32_4, color = 'm')
plt.plot(time_26_4,med_beq, color = 'grey', linestyle = '--',linewidth = 8, label = ('median'))
#plt.plot(time_26_4,mean_beq, color = 'black', linestyle = '-.',linewidth = 8, label = ('mean'))
p = plt.axvspan(19.9, 20, facecolor='green', alpha=0.5, label = ('sunset'))
plt.legend(loc = 2, fontsize = 'medium')
plt.xlabel('Time (Hours of Day)')
plt.ylabel('Becquerel')
plt.ylim(0.2,0.5)
plt.show()