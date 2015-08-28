import matplotlib.pyplot as plt
import numpy as np
import os,csv,sys
fig1 = plt.figure()
fig2 = plt.figure()
from scipy.stats import spearmanr
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

### DETERMINING GEIGER PLATEAU ###
#volt1 = [200,300,400,500,600,700,750,800,850,900,950,1000,1100,800,950]
#count1 = [0,0,0,0,0,422,638,571,694,676,735,831,1077,684,731]
#volt2 = [600,650,700,750,800,850,900,950,1000,1050]
#count2 = [0,0,111,151,511,516,508,511,552,742]
#volt3 = [900,900,850,800,800,750,700,600]
#count3 = [468,372,368,436,308,352,340,0]
#volt4 = [900,900,850,800,800,750,700,600]
#count4 = [508,372,408,408,384,408,268,0]
#plt.plot(volt1,count1,'r^')
#plt.plot(volt2,count2, 'b*')
#plt.plot(volt3,count3, 'cs')
#plt.plot(volt4,count4, 'mo')
#plt.show()

#sys.exit('stop')


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
    
    
### LINEAR SOLVER ###  

def lin_solver(component,beq):
    a = np.vstack([component, np.ones(len(component))]).T
    ans = np.asarray(np.linalg.lstsq(a, beq))
    print "lin fit" ,ans[0]
    return ans  

print 'pressure fit' 
press_ans = lin_solver(pressure,beq_norm)
print 'temperature fit' 
temp_ans = lin_solver(temp,beq_norm)
print 'humidity fit' 
humid_ans = lin_solver(humidity,beq_norm)
print 'other fits'
press_ansa = lin_solver(pressure1,beq1)
press_ansb = lin_solver(pressure2,beq2)
press_ansc = lin_solver(pressure3,beq3)
press_ansd = lin_solver(pressure4,beq4)
temp_ansa = lin_solver(temp1,beq1)
temp_ansb = lin_solver(temp2,beq2)
temp_ansc = lin_solver(temp3,beq3)
temp_ansd = lin_solver(temp4,beq4)
humid_ansa = lin_solver(humidity1,beq1)
humid_ansb = lin_solver(humidity2,beq2)
humid_ansc = lin_solver(humidity3,beq3)
humid_ansd = lin_solver(humidity4,beq4)






#rho_p, P_p = spearmanr(pressure,beq_norm)
#print 'rho pressure:',rho_p
#print 'probability pressure:',round((1-P_p)*100,1),'%'
#rho_t, P_t = spearmanr(temp,beq_norm)
#print 'rho temperature:',rho_t
#print 'probability temperature:',round((1-P_t)*100,1),'%'
#rho_h, P_h = spearmanr(humidity,beq_norm)
#print 'rho humidity:',rho_h
#print 'probability humidity:',round((1-P_h)*100,1),'%'


rho_p, P_p = spearmanr(pressure,beq_norm)
print 'rho pressure:',rho_p
print 'probability pressure:',P_p
rho_t, P_t = spearmanr(temp,beq_norm)
print 'rho temperature:',rho_t
print 'probability temperature:',P_t
rho_h, P_h = spearmanr(humidity,beq_norm)
print 'rho humidity:',rho_h
print 'probability humidity:',P_h




### PLOTTING ###
plt.rc('font', family='serif')
loc_day = np.where(np.logical_and(time <=  19.75, time >=  5))
loc_night = np.where(np.logical_or(time <=  5,time >  19.75))


plt1 = fig1.add_subplot(1,1,1)
d = plt1.plot(time[loc_day],beq_norm[loc_day],'y.',label = ('day'))
n = plt1.plot(time[loc_night],beq_norm[loc_night],'k.',label = ('night')) 
#fig1.legend((d, n), ('day', 'night'), loc = 4)
#plt1.plot(time1,beq1_norm,'r.')
#plt1.plot(time2,beq2_norm,'b.')
#plt1.plot(time3,beq3_norm,'c.')
#plt1.plot(time4,beq4_norm,'m.')
plt1.set_xlabel('Time of Day (Hours)')
plt1.set_ylabel('Radioactivity (Bq)')
plt1.set_xlim((min(time)-.05*min(time)),(max(time)+.05*max(time)))

#for i in range(len(loc_day[0])):
beq_day = beq_norm[loc_day]    
#for i in range(len(loc_night[0])):
beq_night = beq_norm[loc_night]
avg_day = np.mean(beq_day)
avg_night = np.mean(beq_night)
p1 = plt1.axvspan(19.70, 20, facecolor='red', alpha=0.5, label = ('sunrise/sunset'))
p2 = plt1.axvspan(6.1, 6.4, facecolor='red', alpha=0.5)
l1 = plt1.axhline(y=avg_day, color = 'y', label = ('day mean'))
l2 = plt1.axhline(y=avg_night, color = 'k', label = ('night mean'))
plt1.legend(loc = 'upper center', fontsize = 'medium')




plt2 = fig2.add_subplot(3,1,1)
plt2.plot(pressure[loc_day],beq_norm[loc_day],'y.',label = ('day'))
plt2.plot(pressure[loc_night],beq_norm[loc_night],'k.',label = ('night'))
#plt2.plot(pressure1,beq1_norm,'r.')
#plt2.plot(pressure2,beq2_norm,'b.')
#plt2.plot(pressure3,beq3_norm,'c.')
#plt2.plot(pressure4,beq4_norm,'m.')
plt2.set_xlabel('Pressure (hPa)')
plt2.set_ylabel('Radioactivity (Bq)')
plt2.set_xlim((min(pressure)-.0005*min(pressure)),(max(pressure)+.0005*max(pressure)))
press_ans = press_ans[0]
x2 = np.arange((min(pressure)-.0005*min(pressure)),(max(pressure)+.005*max(pressure)))
def y2(x):
    return press_ans[0]*x+press_ans[1]
plt2.plot(x2, y2(x2),label = ('linear fit'))
#plt2.legend(loc = 3, fontsize = 'medium')


plt3 = fig2.add_subplot(3,1,2)
plt3.plot(temp[loc_day],beq_norm[loc_day],'y.',label = ('day'))
plt3.plot(temp[loc_night],beq_norm[loc_night],'k.',label = ('night'))
#plt3.plot(temp1,beq1_norm,'r.')
#plt3.plot(temp2,beq2_norm,'b.')
#plt3.plot(temp3,beq3_norm,'c.')
#plt3.plot(temp4,beq4_norm,'m.')
plt3.set_xlabel('Temperature ($^\circ$C)')
plt3.set_ylabel('Radioactivity (Bq)')
plt3.set_xlim((min(temp)-.05*min(temp)),(max(temp)+.05*max(temp)))
temp_ans = temp_ans[0]
x3 = np.arange((min(temp)-.05*min(temp)),(max(temp)+.15*max(temp)))
def y3(x):
    return temp_ans[0]*x+temp_ans[1]
plt3.plot(x3, y3(x3),label = ('linear fit'))
plt3.legend(loc = 2, fontsize = 'medium')


plt4 = fig2.add_subplot(3,1,3)
plt4.plot(humidity[loc_day],beq_norm[loc_day],'y.',label = ('day'))
plt4.plot(humidity[loc_night],beq_norm[loc_night],'k.',label = ('night'))
#plt4.plot(humidity1,beq1_norm,'r.')
#plt4.plot(humidity2,beq2_norm,'b.')
#plt4.plot(humidity3,beq3_norm,'c.')
#plt4.plot(humidity4,beq4_norm,'m.')
plt4.set_xlabel('Humidity (\%)')
plt4.set_ylabel('Radioactivity (Bq)')
plt4.set_xlim((min(humidity)-.05*min(humidity)),(max(humidity)+.05*max(humidity)))
humid_ans = humid_ans[0]
x4 = np.arange((min(humidity)-.05*min(humidity)),(max(humidity)+.15*max(humidity)))
def y4(x):
    return humid_ans[0]*x+humid_ans[1]
plt4.plot(x4, y4(x4),label = ('linear fit'))
#plt4.legend(loc = 4, fontsize = 'medium')

fig2.text(.4, .11, r'$\rho = 0.13$', fontsize=15)
fig2.text(.4, .395, r'$\rho = 0.09$', fontsize=15)
fig2.text(.4, .68, r'$\rho = -0.18$', fontsize=15)

fig2.text(.6, .110, r'$\mathit{p\mbox{-}value} = 0.0057$', fontsize=15)
fig2.text(.6, .395, r'$\mathit{p\mbox{-}value} = 0.1961$', fontsize=15)
fig2.text(.6, .685, r'$\mathit{p\mbox{-}value} = 0.0422$', fontsize=15)


#
#fig2.legend((d, n), ('day', 'night'), 'upper right')


#total_steradian = 2*np.pi
#sun_steradian = 6.87*(10**(-5)) 
#steradian_fraction = sun_steradian/total_steradian
#print 'steradian fraction:'
#print steradian_fraction
#becquerel_fraction = (avg_night-avg_day)/avg_night
#print 'becquerel fraction:'
#print becquerel_fraction






'''

figa = plt.figure()
figb = plt.figure()
figc = plt.figure()
figd = plt.figure()


plt1a = figa.add_subplot(2,2,1) 
plt1a.plot(time1,beq1,'r^')
#plt1a.plot(time2,beq2_norm,'b*')
#plt1a.plot(time3,beq3_norm,'cs')
#plt1a.plot(time4,beq4_norm,'mo')
plt1a.set_xlabel('time of day (hours)')
plt1a.set_ylabel('becquerel')
plt1a.set_xlim((min(time)-.05*min(time)),(max(time)+.05*max(time)))

plt2a = figa.add_subplot(2,2,2)
plt2a.plot(pressure1,beq1,'r^')
#plt2a.plot(pressure2,beq2_norm,'b*')
#plt2a.plot(pressure3,beq3_norm,'cs')
#plt2a.plot(pressure4,beq4_norm,'mo')
plt2a.set_xlabel('pressure (hPa)')
plt2a.set_ylabel('becquerel')
plt2a.set_xlim((min(pressure)-.0005*min(pressure)),(max(pressure)+.0005*max(pressure)))
press_ans = press_ansa[0]
x2 = np.arange(min(pressure),max(pressure))
def y2a(x):
    return press_ansa[0]*x+press_ansa[1]
plt2a.plot(x2, y2(x2))

plt3a = figa.add_subplot(2,2,3)
plt3a.plot(temp1,beq1,'r^')
#plt3.plot(temp2,beq2_norm,'b*')
#plt3.plot(temp3,beq3_norm,'cs')
#plt3.plot(temp4,beq4_norm,'mo')
plt3a.set_xlabel('temperature (C)')
plt3a.set_ylabel('becquerel')
plt3a.set_xlim((min(temp)-.05*min(temp)),(max(temp)+.05*max(temp)))
temp_ans = temp_ansa[0]
x3 = np.arange(min(temp),max(temp))
def y3a(x):
    return temp_ansa[0]*x+temp_ansa[1]
plt3a.plot(x3, y3(x3))


plt4a = figa.add_subplot(2,2,4)
plt4a.plot(humidity1,beq1,'r^')
#plt4.plot(humidity2,beq2_norm,'b*')
#plt4.plot(humidity3,beq3_norm,'cs')
#plt4.plot(humidity4,beq4_norm,'mo')
plt4a.set_xlabel('humidity (%)')
plt4a.set_ylabel('becquerel')
plt4a.set_xlim((min(humidity)-.05*min(humidity)),(max(humidity)+.05*max(humidity)))
humid_ans = humid_ansa[0]
x4 = np.arange(min(humidity),max(humidity))
def y4a(x):
    return humid_ansa[0]*x+humid_ansa[1]
plt4a.plot(x4, y4(x4))

loc_day = np.where(time1 <  19)
for i in range(len(loc_day[0])):
    beq_day = beq1[loc_day[0:i]]    
loc_night = np.where(time1 >=  19)
for i in range(len(loc_night[0])):
    beq_night = beq1[loc_night[0:i]]
avg_day = np.mean(beq_day)
avg_night = np.mean(beq_night)

p = plt1a.axvspan(19.70, 20, facecolor='g', alpha=0.5)
l1 = plt1a.axhline(y=avg_day, color = 'y')
l2 = plt1a.axhline(y=avg_night, color = 'k')








plt1b = figb.add_subplot(2,2,1) 
plt1b.plot(time2,beq2,'b*')
#plt1a.plot(time3,beq3_norm,'cs')
#plt1a.plot(time4,beq4_norm,'mo')
plt1a.set_xlabel('time of day (hours)')
plt1a.set_ylabel('becquerel')
plt1a.set_xlim((min(time)-.05*min(time)),(max(time)+.05*max(time)))

plt2b = figb.add_subplot(2,2,2)
plt2b.plot(pressure2,beq2,'b*')
#plt2a.plot(pressure3,beq3_norm,'cs')
#plt2a.plot(pressure4,beq4_norm,'mo')
plt2b.set_xlabel('pressure (hPa)')
plt2b.set_ylabel('becquerel')
plt2b.set_xlim((min(pressure)-.0005*min(pressure)),(max(pressure)+.0005*max(pressure)))
press_ans = press_ansb[0]
x2 = np.arange(min(pressure),max(pressure))
def y2b(x):
    return press_ansb[0]*x+press_ansb[1]
plt2b.plot(x2, y2(x2))

plt3b = figb.add_subplot(2,2,3)
plt3b.plot(temp2,beq2,'b*')
#plt3.plot(temp3,beq3_norm,'cs')
#plt3.plot(temp4,beq4_norm,'mo')
plt3b.set_xlabel('temperature (C)')
plt3b.set_ylabel('becquerel')
plt3b.set_xlim((min(temp)-.05*min(temp)),(max(temp)+.05*max(temp)))
temp_ans = temp_ansa[0]
x3 = np.arange(min(temp),max(temp))
def y3b(x):
    return temp_ansb[0]*x+temp_ansb[1]
plt3b.plot(x3, y3(x3))

plt4b = figb.add_subplot(2,2,4)
plt4b.plot(humidity2,beq2,'b*')
#plt4.plot(humidity3,beq3_norm,'cs')
#plt4.plot(humidity4,beq4_norm,'mo')
plt4b.set_xlabel('humidity (%)')
plt4b.set_ylabel('becquerel')
plt4b.set_xlim((min(humidity)-.05*min(humidity)),(max(humidity)+.05*max(humidity)))
humid_ans = humid_ansb[0]
x4 = np.arange(min(humidity),max(humidity))
def y4b(x):
    return humid_ansb[0]*x+humid_ansb[1]
plt4b.plot(x4, y4(x4))

loc_day = np.where(time2 <  19)
for i in range(len(loc_day[0])):
    beq_day = beq2[loc_day[0:i]]    
loc_night = np.where(time2 >=  19)
for i in range(len(loc_night[0])):
    beq_night = beq2[loc_night[0:i]]
avg_day = np.mean(beq_day)
avg_night = np.mean(beq_night)

p = plt1b.axvspan(19.70, 20, facecolor='g', alpha=0.5)
l1 = plt1b.axhline(y=avg_day, color = 'y')
l2 = plt1b.axhline(y=avg_night, color = 'k')







plt1c = figc.add_subplot(2,2,1) 
plt1c.plot(time3,beq3,'cs')
#plt1a.plot(time4,beq4_norm,'mo')
plt1c.set_xlabel('time of day (hours)')
plt1c.set_ylabel('becquerel')
plt1c.set_xlim((min(time)-.05*min(time)),(max(time)+.05*max(time)))

plt2c = figc.add_subplot(2,2,2)
plt2c.plot(pressure3,beq3,'cs')
#plt2a.plot(pressure4,beq4_norm,'mo')
plt2c.set_xlabel('pressure (hPa)')
plt2c.set_ylabel('becquerel')
plt2c.set_xlim((min(pressure)-.0005*min(pressure)),(max(pressure)+.0005*max(pressure)))
press_ans = press_ansc[0]
x2 = np.arange(min(pressure),max(pressure))
def y2c(x):
    return press_ansc[0]*x+press_ansc[1]
plt2c.plot(x2, y2(x2))

plt3c = figc.add_subplot(2,2,3)
plt3c.plot(temp3,beq3,'cs')
#plt3.plot(temp4,beq4_norm,'mo')
plt3c.set_xlabel('temperature (C)')
plt3c.set_ylabel('becquerel')
plt3c.set_xlim((min(temp)-.05*min(temp)),(max(temp)+.05*max(temp)))
temp_ans = temp_ansc[0]
x3 = np.arange(min(temp),max(temp))
def y3c(x):
    return temp_ansc[0]*x+temp_ansc[1]
plt3c.plot(x3, y3(x3))

plt4c = figc.add_subplot(2,2,4)
plt4c.plot(humidity3,beq3,'cs')
#plt4.plot(humidity4,beq4_norm,'mo')
plt4c.set_xlabel('humidity (%)')
plt4c.set_ylabel('becquerel')
plt4c.set_xlim((min(humidity)-.05*min(humidity)),(max(humidity)+.05*max(humidity)))
humid_ans = humid_ansc[0]
x4 = np.arange(min(humidity),max(humidity))
def y4c(x):
    return humid_ansc[0]*x+humid_ansc[1]
plt4c.plot(x4, y4(x4))

loc_day = np.where(time3 <  19)
for i in range(len(loc_day[0])):
    beq_day = beq3[loc_day[0:i]]    
loc_night = np.where(time3 >=  19)
for i in range(len(loc_night[0])):
    beq_night = beq3[loc_night[0:i]]
avg_day = np.mean(beq_day)
avg_night = np.mean(beq_night)

p = plt1c.axvspan(19.70, 20, facecolor='g', alpha=0.5)
l1 = plt1c.axhline(y=avg_day, color = 'y')
l2 = plt1c.axhline(y=avg_night, color = 'k')







plt1d = figd.add_subplot(2,2,1) 
plt1d.plot(time4,beq4,'mo')
plt1d.set_xlabel('time of day (hours)')
plt1d.set_ylabel('becquerel')
plt1d.set_xlim((min(time)-.05*min(time)),(max(time)+.05*max(time)))

plt2d = figd.add_subplot(2,2,2)
plt2d.plot(pressure4,beq4,'mo')
plt2d.set_xlabel('pressure (hPa)')
plt2d.set_ylabel('becquerel')
plt2d.set_xlim((min(pressure)-.0005*min(pressure)),(max(pressure)+.0005*max(pressure)))
press_ans = press_ansd[0]
x2 = np.arange(min(pressure),max(pressure))
def y2d(x):
    return press_ansd[0]*x+press_ansd[1]
plt2d.plot(x2, y2(x2))

plt3d = figd.add_subplot(2,2,3)
plt3d.plot(temp4,beq4,'mo')
plt3d.set_xlabel('temperature (C)')
plt3d.set_ylabel('becquerel')
plt3d.set_xlim((min(temp)-.05*min(temp)),(max(temp)+.05*max(temp)))
temp_ans = temp_ansd[0]
x3 = np.arange(min(temp),max(temp))
def y3d(x):
    return temp_ansd[0]*x+temp_ansd[1]
plt3d.plot(x3, y3(x3))

plt4d = figd.add_subplot(2,2,4)
plt4d.plot(humidity4,beq4,'mo')
plt4d.set_xlabel('humidity (%)')
plt4d.set_ylabel('becquerel')
plt4d.set_xlim((min(humidity)-.05*min(humidity)),(max(humidity)+.05*max(humidity)))
humid_ans = humid_ansd[0]
x4 = np.arange(min(humidity),max(humidity))
def y4d(x):
    return humid_ansd[0]*x+humid_ansd[1]
plt4d.plot(x4, y4(x4))

loc_day = np.where(time4 <  19)
for i in range(len(loc_day[0])):
    beq_day = beq4[loc_day[0:i]]    
loc_night = np.where(time4 >=  19)
for i in range(len(loc_night[0])):
    beq_night = beq4[loc_night[0:i]]
avg_day = np.mean(beq_day)
avg_night = np.mean(beq_night)

p = plt1d.axvspan(19.70, 20, facecolor='g', alpha=0.5)
l1 = plt1d.axhline(y=avg_day, color = 'y')
l2 = plt1d.axhline(y=avg_night, color = 'k')
'''
plt.show()
