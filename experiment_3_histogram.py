import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from numpy.random import randn
import numpy as np
import os,csv,sys
import math
import matplotlib.patches as mpatches
from matplotlib import rc
from scipy.stats import ttest_ind
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


fig = plt.figure()
plt1 = fig.add_subplot(1,2,1) 
plt2 = fig.add_subplot(1,2,2) 



### READ DATA ###
os.chdir('/Users/Thomas/Desktop/Code/Python/School/Diurnal_Cosmic_Radiation')
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

loc_day = np.where(np.logical_and(time <=  19.75, time >=  5))
loc_night = np.where(np.logical_or(time <=  5,time >  19.75))
#for i in range(len(loc_day[0])):
beq_day = beq_norm[loc_day]    
#for i in range(len(loc_night[0])):
beq_night = beq_norm[loc_night]

hist_night, bins_night = np.histogram(beq_night, bins = 25)
width_night = 0.7 * (bins_night[1] - bins_night[0])
center_night = (bins_night[:-1] + bins_night[1:]) / 2
hist_day, bins_day = np.histogram(beq_day, bins = 25)
width_day = 0.7 * (bins_day[1] - bins_day[0])
center_day = (bins_day[:-1] + bins_day[1:]) / 2
hist, bins = np.histogram(beq, bins = 25)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

mean_n1 = np.mean(beq_night)
variance_n = np.var(beq_night)
sigma_n1 = math.sqrt(variance_n)

#x = np.linspace(0.2,0.5,100)
#plt1.plot(x,mlab.normpdf(x,mean_n,sigma_n), color = 'b',linewidth = 2)

mean_d1 = np.mean(beq_day)
variance_d = np.var(beq_day)
sigma_d1 = math.sqrt(variance_d)

#x = np.linspace(0.2,0.5,100)
#plt1.plot(x,mlab.normpdf(x,mean_d,sigma_d), color = 'y',linewidth = 2)
plt.rc('font', family='serif')
plt1.bar(center-width/2, hist_day, align='center', width=width/2 ,color = 'yellow', label = ('day; mean:'+str(round(mean_d1,3))))
plt1.bar(center, hist_night, align='center', width=width/2, color = 'cyan', label = ('night; mean:'+str(round(mean_n1,3))))
#p = plt1.axvspan(mean_n-0.003, mean_n+0.003, facecolor='blue', alpha=0.5)
#p = plt1.axvspan(mean_d-0.003, mean_d+0.003, facecolor='y', alpha=0.5)
plt1.legend(loc = 2, fontsize = 'medium')
plt1.set_xlabel('Observed Radioactivity (Bq)')
plt1.set_ylabel('Frequency')

count = 0
diff_mean = mean_n1-mean_d1
diff_rand = 0
while (diff_rand) <= (mean_n1-mean_d1):
#for i in range(0,999999):
    mu = (mean_d1+mean_n1)/2
    sigma = (sigma_d1+ sigma_n1)/2
    dist_night = sigma * randn(1, len(beq_night)) + mu
    dist_day = sigma * randn(1, len(beq_day)) + mu
    
    mean_d = np.mean(dist_day)
    mean_n = np.mean(dist_night)
    variance_n = np.var(dist_day)
    sigma_n = math.sqrt(variance_n)
    variance_d = np.var(dist_day)
    sigma_d = math.sqrt(variance_d)
    
    if mean_d >= mean_n:
        diff_rand = mean_d-mean_n
    if mean_n > mean_d:
        diff_rand = mean_n-mean_d
    if diff_rand > diff_mean:
        count +=1
print count

hist_night2, bins_night2 = np.histogram(dist_night, bins = 25)
width_night2 = 0.7 * (bins_night2[1] - bins_night2[0])
center_night2 = (bins_night2[:-1] + bins_night2[1:]) / 2
hist_day2, bins_day2 = np.histogram(dist_day, bins = 25)
width_day2 = 0.7 * (bins_day2[1] - bins_day2[0])
center_day2 = (bins_day2[:-1] + bins_day2[1:]) / 2
hist2, bins2 = np.histogram(beq, bins = 25)
width2 = 0.7 * (bins2[1] - bins2[0])
center2 = (bins2[:-1] + bins2[1:]) / 2

#x2n = np.linspace(0.2,0.5,100)
#plt2.plot(x,mlab.normpdf(x2n,mean_n,sigma_n), color = 'b',linewidth = 2)
#x2d = np.linspace(0.2,0.5,100)
#plt2.plot(x,mlab.normpdf(x2d,mean_d,sigma_d), color = 'y',linewidth = 2)
plt2.bar(center-width/2, hist_day2, align='center', width=width/2 ,color = 'yellow', label = ('day; mean:'+str(round(mean_d,3))))
plt2.bar(center, hist_night2, align='center', width=width/2, color = 'cyan', label = ('night; mean:'+str(round(mean_n,3))))
#p = plt2.axvspan(mean_n-0.003, mean_n+0.003, facecolor='blue', alpha=0.5)
#p = plt2.axvspan(mean_d-0.003, mean_d+0.003, facecolor='y', alpha=0.5)
plt2.legend(loc = 2, fontsize = 'medium')
plt2.set_xlabel('Randomly Generated Radioactivity (Bq)')
plt2.set_ylabel('Frequency')
plt2.set_ylim(0,14)

plt.show()




#t = (np.mean(beq_night)-np.mean(beq_day))/(math.sqrt((((np.std(beq_night))**2)/len(beq_night))+(((np.std(beq_day))**2)/len(beq_day))))
t = ttest_ind(beq_night, beq_day)
print 't-test result:', t