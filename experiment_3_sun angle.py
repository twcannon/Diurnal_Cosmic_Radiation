import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import os,csv,sys
import math
from numpy import sin, cos, arccos 
from scipy.stats import spearmanr


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
date = np.concatenate((date1,date2,date3,date4))


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

    
    
chas_lat = 32.7833
dec = 65.45
freq = (np.pi/24)*2
amp =  dec
shift = (np.pi/24)-12.25
drop = chas_lat
x = np.arange(0,24)

def y(x,amp,freq,shift,drop):
    return amp*sin((x*freq)+shift)+drop
    
def lin_solver(component,beq):
    a = np.vstack([component, np.ones(len(component))]).T
    ans = np.asarray(np.linalg.lstsq(a, beq))
    print "lin fit" ,ans[0]
    return ans  

        
declination = []
for i in range(0, len(time)):
    val = y(time[i],amp,freq,shift,drop)
    declination.append(val)

#def lin_solver(component,beq):
#    a = np.vstack([component, np.ones(len(component))]).T
#    ans = np.asarray(np.linalg.lstsq(a, beq))
#    print "lin fit" ,ans[0]
#    return ans  
#ans = lin_solver(declination,beq_norm)
#ans = ans[0]
#x = np.arange(min(declination),max(declination))
#def y(x):
#    return ans[0]*x+ans[1]
#plt.plot(x, y(x),label = ('linear fit'))
#plt.legend(loc = 2, fontsize = 'small')


rho, P = spearmanr(beq_norm,declination)
print 'rho:',rho
print 'probability:',round((1-P)*100,1),'%'


plt.plot(declination,beq_norm,'g.') 
plt.ylabel('Becquerel')
plt.xlabel('Sun Angle Above Horizon(degrees)')
ans = lin_solver(declination,beq_norm)
ans = ans[0]
x2 = np.arange(-40,100)
def y2(x2):
    return ans[0]*x2+ans[1]
plt.plot(x2, y2(x2),label = ('linear fit'))
#plt.xlim(0,110)
plt.show()