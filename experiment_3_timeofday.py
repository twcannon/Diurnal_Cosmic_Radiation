import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import os,csv,sys
import math
from numpy import sin
from scipy.stats import spearmanr
from scipy.signal import detrend


fig = plt.figure()


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
time = np.concatenate((time1,time2,time3,time4))
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

unique_time = []
for i in range(1,len(time)):
    if time[i] != time[i-1]:
        unique_time.append(round(time[i-1],2))
unique_time.append(round(time[len(time)-1],2))

def f1(seq):
   # not order preserving
   set = {}
   map(set.__setitem__, seq, [])
   return set.keys()
   
unique_time = f1(unique_time)

mean_array = []
std_arr = []
for i in range(0,len(unique_time)):
    loc = np.where(time == unique_time[i])
    mean_beq = np.mean(beq_norm[loc])
    std_beq = (np.std(beq_norm[loc]))
    mean_array.append(mean_beq)
    std_arr.append(std_beq)
for i in range(0,len(unique_time)):   
    if unique_time[i] > 24.:
        unique_time[i] = unique_time[i]-24
        
chas_lat = 32.7833
dec = 43.2
freq = (np.pi/24)*2
amp =  dec
shift = (np.pi/24)-11
drop = chas_lat
x = np.arange(0,30)

def y(x,amp,freq,shift,drop):
    return amp*sin((x*freq)+shift)+drop
    
def lin_solver(component,beq):
    a = np.vstack([component, np.ones(len(component))]).T
    ans = np.asarray(np.linalg.lstsq(a, beq))
    print "lin fit" ,ans[0]
    return ans  

        
declination = []
for i in range(0, len(unique_time)):
    val = y(unique_time[i],amp,freq,shift,drop)
    declination.append(val)
declination_all = []
for i in range(0, len(time)):
    val_all = y(time[i],amp,freq,shift,drop)
    declination_all.append(val_all)
    
    
    
    
humid_ans = lin_solver(humidity,beq_norm)
humid_ans = humid_ans[0]
sub_beq = []
for i in range(0,len(beq_norm)):
    val = beq_norm[i] - ((humidity[i]*humid_ans[0])+humid_ans[1]) 
    sub_beq.append(val)
    
press_ans = lin_solver(pressure,sub_beq)
press_ans = press_ans[0]
sub_sub_beq = []
for i in range(0,len(sub_beq)):
    val = sub_beq[i] - ((pressure[i]*press_ans[0])+press_ans[1]) 
    sub_sub_beq.append(val)
press_new_ans = lin_solver(pressure,sub_sub_beq)


    
        
        
        
#plt1 = fig.add_subplot(3,1,1)
plt2 = fig.add_subplot(2,1,1)
plt3 = fig.add_subplot(2,1,2)
#plt4 = fig.add_subplot(4,1,2)



#plt1.plot(unique_time,mean_array, 'b.')
##plt1.plot(pressure,beq_norm, 'b.')
##plt1.plot(time,beq_norm, 'b.')
##plt1.errorbar(unique_time, mean_array, xerr=std_arr, fmt='b.')#, ecolor=None, elinewidth=None, capsize=3)
#plt1.set_xlabel('Time of Day (hours)')
#plt1.set_ylabel('Becquerel')
##plt1.set_xlim(0,24)

#plt4.plot(unique_time,mean_array, 'b.')
##plt4.plot(pressure,sub_beq, 'b.')
##plt4.plot(time,beq_norm, 'b.')
##plt4.errorbar(unique_time, mean_array, xerr=std_arr, fmt='b.')#, ecolor=None, elinewidth=None, capsize=3)
#plt4.set_xlabel('Time of Day (hours)')
#plt4.set_ylabel('Becquerel')
##plt4.set_xlim(0,24)

plt2.set_xlim(0,24)
plt2.set_xlabel('Time of Day (Hours)')
plt2.set_ylabel('Declination (Degrees)')
plt2.plot(x, y(x,amp,freq,shift,drop),label = ('Height of Cluster Center'),color = 'r', linewidth = 5)
plt2.plot(x, y(x,amp,freq,shift,drop-17),label = ('Height of Cluster Center'),color = 'b')
plt2.plot(x, y(x,amp,freq,shift,drop+17),label = ('Height of Cluster Center'),color = 'b')
l = plt2.axhspan(-40, 0, facecolor='k', alpha=0.5)
plt2.set_ylim(-40,90)


rho, P = spearmanr(declination_all,beq_norm)
#rho, P = spearmanr(declination,mean_array)
#rho, P = spearmanr(declination_all,sub_sub_beq)
print 'rho:',rho
print 'probability:',round((1-P)*100,1),'%'
print 'probability:',P
#plt3.plot(declination,mean_array,'g.')
#ans = lin_solver(declination,mean_array)
#plt3.plot(declination_all,beq_norm,'g.')
#ans = lin_solver(declination_all,beq_norm)
plt3.plot(declination_all,beq_norm,'g.')
ans = lin_solver(declination_all,beq_norm)
ans = ans[0]
x2 = np.arange(-20,100)
def y2(x2):
    return ans[0]*x2+ans[1]
plt3.set_xlim(0,90)
plt3.plot(x2, y2(x2),label = ('linear fit'))
plt3.set_xlabel('Declination of Hotspot (Degrees)')
plt3.set_ylabel('Radioactivity (Bq)')

fig.text(.3, .110, r'$\rho = 0.15$', fontsize=15)
fig.text(.6, .110, r'$\mathit{p\mbox{-}value} = 0.0231$', fontsize=15)
fig.text(.15, .65, r'Horizon', fontsize=15)
fig.text(.4, .8, r'Path of Hotspot', fontsize=15)

plt.show()
