# code copied from AFM3.python on 02/01/17

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 10:26:02 2017

@author: sp4009
"""
'''
AFManalysis 3 imports data as text file

AFM images should be same size e.g. 5um x 5um
exported images should be saved at PNG 300 dpi, 10 inch x 10 inch, NO extras
threshold images formated in NanoScope Analysis with threshold height saved in txt file in same folder
crosssection taken in NanoScope Analysis at any random point
full roughness analyssi exported from Nanoscope Anlaysis
(information lost in image alone so if threshold / crossection / histogram analsysis done using python  results are not accurate)

AFM object created with all aspects. sections below then used to plot each.
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import fftpack
import glob
import os, errno
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
'''
create AFM class
'''
class AFM(object):
    def __init__(self, name, path, scanSizeX = 256, scanSizeY = 256, size_um = 5):
        self.data = pd.read_csv(path, encoding = 'ANSI')
        self.data.columns = [name]
        self.matrix = np.array([self.data[name].values.tolist()[i:i + scanSizeX] for i in range(0, len(self.data), scanSizeY)])
        self.name = name
        self.concentration = int(name[1:3])
        '''
        histogram, cross section, roughness, colour, condition
        '''
        self.hist = np.histogram(self.data, bins = np.linspace(-10, 10, 100), density = True)[0]
        self.cross = self.matrix[scanSizeY/2]
        self.roughRMS = np.sqrt(np.mean(np.square(self.data.values)))         
        self.colour, self.condition = self.f_colourAndCondition()
        self.fourier, self.fourierProfile = self.f_four(self.matrix)
        self.thresh, self.threshRatio = self.threshold(self.matrix)
        
        self.matrixc = self.correction()
        self.histc = np.histogram(self.matrixc, bins = np.linspace(-10, 10, 100), density = True)[0]
        self.crossc = self.matrixc[scanSizeY/2]
        self.fourierc, self.fourierProfilec = self.f_four(self.matrixc)
        self.threshc, self.threshRatioc = self.threshold(self.matrixc)
        
    def correction(self):
        maxHeight = np.linspace(-10, 10, 100)[np.where(self.hist == max(self.hist))] # several max heights may be present so take first
        return self.matrix - maxHeight[0]
         
    
        
    def f_colourAndCondition(self):
        if self.name[0] == 'A':
            c = 'k-'
            con = 'dark - in-situ'
        elif self.name[0] == 'B':
            c = 'k--'
            con = 'dark - steady state'
        elif self.name[0] == 'C':
            c = 'r-'
            con = 'light - in-situ'
        elif self.name[0] == 'D':
            c = 'r--'
            con = 'light - steady state'
        else:
            c = 'b'
            con = 'unkown'
        return c, con        

    def threshold(self, matrix, height = 2.5):
        thresh = matrix.copy()
        thresh[thresh <= height] = 0
        thresh[thresh > height] = 255
        return thresh, thresh[thresh>0].size/thresh.size

    def f_four(self, matrix):
        '''
        returns fourier transform of imported image and calculates the radial profile
        '''
        F1 = fftpack.fft2(matrix) # create fft
        F2 = abs(fftpack.fftshift(F1)) # shift fft in F1 to centre and abs values
        centre = [F2.shape[0] / 2, F2.shape[1] / 2] #find centre of image
        x, y = np.indices((F2.shape)) #return array of indicies [0,0..][1,1...] and [1,2,3...][1,2,3..]
        r = np.sqrt((x - centre[0])**2 + (y - centre[1])**2).astype(np.int) # get distances of x, y, from centre
        radialprofile = np.bincount(r.ravel(), F2.ravel()) / np.bincount(r.ravel()) # sum F2 at distances / distances   
        return F2, radialprofile
    
        
        
        #self.profT = pd.read_excel(path + 'ThickAndRoughData.xlsx')[self.name.lower()][3]          #profilometerThickness
        #self.profPa = pd.read_excel(path + 'ThickAndRoughData.xlsx')[self.name.lower()][10] 
        


#%%
'''
120degC heating - load data
'''
AFM120D80 = {}
timesD80 = []
path = 'H:/Research/OPV/Neutron Scattering/AFM/120degCheating/D80/cropped/txt/'

for path in np.sort(glob.glob(path +  "*.txt")):
    name = path[70:-24]
    AFM120D80[name] = AFM(name, path)
    timesD80.append(int(name[4:-3]))

AFM120B80 = {}
timesB80 = []
path = 'H:/Research/OPV/Neutron Scattering/AFM/120degCheating/B80/cropped/txt/'

for path in np.sort(glob.glob(path +  "*.txt")):
    name = path[70:-24]
    AFM120B80[name] = AFM(name, path)
    timesB80.append(int(name[4:-3]))    

Bkeys = [k for k in AFM120B80.keys()]
Dkeys = [k for k in AFM120D80.keys()]
'''    
manual corrections
'''
AFM120B80['B80_0115min'].matrixc = AFM120B80['B80_0115min'].matrixc - 1.5
    
#%% creating dataframes
fourD80 = pd.concat(((pd.DataFrame(AFM120D80[i].fourierProfile, columns=[int(i[4:-3])])) for i in sorted(AFM120D80.keys())), axis=1)
fourD80.index = np.linspace(0,5*len(fourD80)/256, len(fourD80))
histD80 = pd.concat(((pd.DataFrame(AFM120D80[i].histc, columns=[int(i[4:-3])])) for i in sorted(AFM120D80.keys())), axis=1)
histD80.index = np.linspace(-10,10, len(histD80))
crossD80 = pd.concat(((pd.DataFrame(AFM120D80[i].cross, columns=[int(i[4:-3])])) for i in sorted(AFM120D80.keys())), axis=1)
crossD80.index = np.linspace(0,5, len(crossD80))

fourB80 = pd.concat(((pd.DataFrame(AFM120B80[i].fourierProfile, columns=[int(i[4:-3])])) for i in sorted(AFM120B80.keys())), axis=1)
fourB80.index = np.linspace(0,5*len(fourB80)/256, len(fourB80))
histB80 = pd.concat(((pd.DataFrame(AFM120B80[i].histc, columns=[int(i[4:-3])])) for i in sorted(AFM120B80.keys())), axis=1)
histB80.index = np.linspace(-10,10, len(histB80))
crossB80 = pd.concat(((pd.DataFrame(AFM120B80[i].cross, columns=[int(i[4:-3])])) for i in sorted(AFM120B80.keys())), axis=1)
crossB80.index = np.linspace(0,5, len(crossB80))

#%% plotting AFMs
exp = AFM120B80
plt.figure(figsize = (16,2))
gs1 = gridspec.GridSpec(1, 8)
gs1.update(wspace=0, hspace=0) 
j = 0
for i in Bkeys[:5]:
    ax1 = plt.subplot(gs1[j])
    ax1.imshow(exp[i].matrixc, 'afmhot', vmin = -10, vmax = 10)
    ax1.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
    ax1.set_title(str(int(i[4:-3])) + ' min')
    ax1.axis('off')
    ax2 = plt.subplot(gs1[j])
    ax2 = inset_axes(ax1, width="40%", height='40%', loc=3)
    ax2.imshow(exp[i].fourier, 'afmhot', vmin = 0, vmax = 10000)
    ax2.set_xlim(100,156)
    ax2.set_ylim(100,156)
    ax2.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
    ax2.axis('off')
    j += 1


#%% plotting fourier profiles
f, (ax1, ax2) = plt.subplots(2,1, figsize = (3,3))
fourB80[fourB80.keys()[[0,1,2,3,4]]].plot(ax = ax1, colormap = 'coolwarm', xlim=(0.03,0.5), ylim=(0,1.2e4), marker='o', mfc='white', ms = 4,  mew=1, lw = 0, legend = False)
fourD80[fourD80.keys()[[0,1,2,3,4]]].plot(ax = ax2, colormap = 'coolwarm', xlim=(0.03,0.5), ylim=(0,1.2e4), marker='o', mfc='white', ms = 4,  mew=1, lw = 0, legend = False)
ax1.set_xticklabels(())
ax1.yaxis.set_ticks([])
ax2.yaxis.set_ticks([])
ax2.set_ylabel('                            2D  Fourier Transform')
ax2.set_xlabel('q (\u03BCm$^{-1}$)')
#ax2.legend(bbox_to_anchor=(1.05, 0.3), loc=4,borderaxespad=1, title = 't (min)', framealpha = 1)


plt.subplots_adjust(hspace = 0)

#%% plotting histogram
f, (ax1, ax2) = plt.subplots(2,1, figsize = (3,3))

histB80[fourB80.keys()[:5]].plot(ax = ax1, colormap = 'coolwarm', xlim=(-10,10), ylim=(0,0.29), marker='o', mfc='white', ms = 3,  mew=1, lw = 1, legend = False)
histD80[fourD80.keys()[:5]].plot(ax = ax2, colormap = 'coolwarm', xlim=(-10,10), ylim=(0,0.29), marker='o', mfc='white', ms = 3,  mew=1, lw = 1, legend = False)
ax1.set_xticklabels(())
#ax1.set_ylabel('Frequency')
ax2.set_ylabel('                         Frequency')
ax2.set_xlabel('Height (nm)')
ax2.legend(bbox_to_anchor=(1.15, -0.1), loc=4,borderaxespad=1, title = 't (min)', framealpha = 1)
#histD80[fourD80.keys()[4]].plot(ax = ax2, color = [0.75,0.157,0.27], xlim=(-10,10), ylim=(0,0.29), marker='o', mfc='white', ms = 3,  mew=2, lw = 1, legend = False)
#histB80[fourB80.keys()[2]].plot(ax = ax1, color = [0.89,0.886,0.886], xlim=(-10,10), ylim=(0,0.29), marker='o', mfc='white', ms = 3,  mew=2, lw = 1, legend = False)

plt.subplots_adjust(hspace = 0)


#%% plotting roughness
Brough = [AFM120B80[k].roughRMS for k in AFM120B80.keys()]
Drough = [AFM120D80[k].roughRMS for k in AFM120D80.keys()]

f, ax = plt.subplots(figsize = (3,3))
ax.plot(timesB80[:6], Brough[:6], 'ko-')
ax.plot(timesD80[:6], Drough[:6], 'ro-')
ax.set_xlabel('annealing time (min)')
ax.set_ylabel('Roughness, Rq (nm)')
#ax.set_xscale('log')
#ax.legend(['Dark', 'Light'])

#%% fourier plots
exp = AFM120D80

f = plt.figure(figsize = (2,2))
plt.imshow(exp[Dkeys[3]].matrixc, 'afmhot', vmin = -10, vmax = 10)
plt.xlim(100,156)
plt.ylim(100,156)
plt.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
plt.colorbar(label='Height (nm)')



#%%

'''
########################################################################################################
initial exp on AFManalysis 3. loading 5 um samples
'''
experiments = ['A', 'B', 'C', 'D']
ratios = ['20', '35', '50', '65', '80']

dic = {}

path = 'H:/Research/OPV/Neutron Scattering/AFM/BestTopographys/txt/'
allFiles = np.sort(glob.glob(path +  "*.txt"))
for path in allFiles:
    name = path[59:62]
    try:
        dic[name] = AFM(name, path)
    except:
        print(name)
        pass

  
#dickeys = [k for k in dic.keys()]
'''
creating dataframes. check index
'''
four = pd.concat(((pd.DataFrame(dic[i].fourierProfile, columns=[i])) for i in sorted(dic.keys())), axis=1)
four.index = np.linspace(0,5*len(fourD80)/256, len(fourD80))
hist = pd.concat(((pd.DataFrame(dic[i].histc, columns=[i])) for i in sorted(dic.keys())), axis=1)
hist.index = np.linspace(-10,10, len(histD80))
cross = pd.concat(((pd.DataFrame(dic[i].cross, columns=[i])) for i in sorted(dic.keys())), axis=1)
cross.index = np.linspace(0,5, len(crossD80))

def exp(exp):
    a = []
    for i in ratios:
        a.append(exp + i)
    return a

#%% plotting AFMs

plt.figure(figsize = (16,2))
gs1 = gridspec.GridSpec(1, 8)
gs1.update(wspace=0, hspace=0) 
j = 0
for i in exp('A'):
    ax1 = plt.subplot(gs1[j])
    ax1.imshow(dic[i].thresh, 'afmhot', vmin = -10, vmax = 10)
    ax1.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
    ax1.set_title(i[1:] + '%')
    ax1.axis('off')
    ax2 = plt.subplot(gs1[j])
    ax2 = inset_axes(ax1, width="40%", height='40%', loc=3)
    ax2.imshow(dic[i].fourier, 'afmhot', vmin = 0, vmax = 10000)
    ax2.set_xlim(100,156)
    ax2.set_ylim(100,156)
    ax2.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
    ax2.axis('off')
    j += 1

#%% plotting histogram
plt.figure(figsize = (16,1))
gs1 = gridspec.GridSpec(1, 8)
gs1.update(wspace=0, hspace=0) 
j = 0
for i in ratios:
    ax = plt.subplot(gs1[j])
    hist['B' + i].plot.area(ax = ax, xlim=(-10,10), color = 'k')#, marker='o', mfc='white', ms = 4,  mew=1, lw = 1, legend = False)
    hist['D' + i].plot.area(ax = ax, xlim=(-10,10), color = 'r')#, marker='o', mfc='white', ms = 4,  mew=1, lw = 1, legend = False)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([-5,0,5])
    #ax.set_yscale('log')
    if j == 0:
        ax.set_ylabel('Frequency')
    if j == 2:
        ax.set_xlabel('Height (nm)') 
    j += 1

#%% plotting fourier
plt.figure(figsize = (16,1))
gs1 = gridspec.GridSpec(1, 8)
gs1.update(wspace=0, hspace=0) 
j = 0
for i in ratios:
    ax = plt.subplot(gs1[j])
    four['A' + i].plot(ax = ax, colormap = 'coolwarm', xlim=(0,1), ylim=(0, 1e4), marker='o', color = 'k', mfc='white', ms = 4,  mew=1, lw = 1, legend = False)
    four['C' + i].plot(ax = ax, colormap = 'coolwarm', xlim=(0,1), ylim=(0, 1e4),  marker='o', color = 'r', mfc='white', ms = 4,  mew=1, lw = 1, legend = False)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([0,0.5])
    if j == 0:
        ax.set_ylabel('2D FFT')
    if j == 2:
        ax.set_xlabel('q (\u03BCm$^{-1}$)') 
    j += 1
#%% plotting roughness

Rough = [dic[k].roughRMS for k in dic.keys()]
'''
Rough values taken from orginial AFM software
'''
Rough[14] = 4.04
Rough[19] = 3.74 # D80

f, ax = plt.subplots(figsize = (3,3))
ax.plot(ratios, Rough[15:], 'rx--')
ax.plot(ratios, Rough[5:10], 'kx--')
ax.plot(ratios, Rough[10:15], 'ro-')
ax.plot(ratios, Rough[:5], 'ko-')

ax.set_xlabel('PCBM conc (%)')
ax.set_ylabel('Roughness, Rq (nm)')

ax.legend(['85\u00B0C Light', '85\u00B0C Dark', '100\u00B0C Light', '100\u00B0C Dark'])

#%%


# concat object hist and cross dataframs
hist = pd.DataFrame([])
cross = pd.DataFrame([])
for i in d.keys():
    hist = pd.concat([hist, d[i].hist], axis = 1)
    cross = pd.concat([cross, d[i].cross], axis = 1)   

#%%  RMS value from Nanoscope and profilometer
f1, ax1 = plt.subplots(figsize=(3,3))
f2, ax2 = plt.subplots(figsize=(3,3))
f3, ax3 = plt.subplots(figsize=(3,3))
f3, ax4 = plt.subplots(figsize=(3,3))
for i in ['B','D']:
    AFMrq = []
    AFMra = []
    profPa = []
    profT = []
    for j in ratios:
        AFMrq.append(d[str(i.lower()+j)].AFMrq)
        AFMra.append(d[str(i.lower()+j)].AFMra)
        profPa.append(d[str(i.lower()+j)].profPa)
        profT.append(d[str(i.lower()+j)].profT)
        
    ax1.set_xlabel('PCBM concentration (%)')
    ax1.plot(ratios, AFMrq, (d[str(i.lower()+j)].colour + 'o'), markerfacecolor='white')
    ax1.set_ylabel('AFM Rq (nm)')
    ax2.set_xlabel('PCBM concentration (%)')
    ax2.plot(ratios, AFMra, (d[str(i.lower()+j)].colour + 'o'), markerfacecolor='white') 
    ax2.set_ylabel('AFM Ra (nm)')
    ax3.set_xlabel('PCBM concentration (%)')    
    ax3.plot(ratios, profPa, (d[str(i.lower()+j)].colour + 'o'), markerfacecolor='white') 
    ax3.set_ylabel('profilometer Pa (nm)')
    ax4.set_xlabel('PCBM concentration (%)')    
    ax4.plot(ratios, profT, (d[str(i.lower()+j)].colour + 'o'), markerfacecolor='white', label = d[str(i.lower()+j)].condition) 
    ax4.set_ylabel('thickness (nm)')
    #ax4.legend(loc=2,bbox_to_anchor=(1, 1))


#%% plotting histrograms
hist.interpolate().plot.area(subplots=True,layout=(4,5), figsize=(10,4), style = ['k-','k-','k-','k-','k-','k-','k-','k-','k-','k-', 'r','r','r','r','r','r','r','r','r','r',], xlim=(-10,10), ylim=(0,2), xticks=[-10,0,10], yticks=[], legend = False, sharex = True, sharey = True)
#%% plotting cross sections
crossEdit = cross + 100
cross.interpolate().plot(subplots=True,layout=(4,5), figsize=(10,4), style = ['k-','k-','k-','k-','k-','k-','k-','k-','k-','k-', 'r','r','r','r','r','r','r','r','r','r',], xlim=(0,5), ylim=(-10,10), yticks=[-10,0,10], xticks=[], legend = False, sharex = True, sharey = True)
#%% plot ratio of bw images and from height of scale in nanoscope analysis with histogram. turns out image analysis not good.
#threshold
f = plt.figure(figsize=(3,3))
for i in experiments:
    x = [] # list of ratios, 'j'
    y = [] # 
    z = []
    for j in ratios:
        if not d[str(i.lower()+j)].imbw == None:
            t = np.array(np.asarray(d[str(i.lower()+j)].imbw)) # convert bwim to an array
            q = (d[str(i.lower()+j)].hist[:d[str(i.lower()+j)].imbwHeight].sum() / d[str(i.lower()+j)].hist.sum())[0]  # sum array < value used to make threshold images in naaoscope. This therefore does not assume even distribution
            x.append(j)
            y.append(t[t < 2].size/t.size) # select total size of array black (pixel < 0) over total array size
            z.append(q)
    #plt.plot(x, y, (d[str(i.lower()+j)].colour + 'x'))
    plt.plot(x, z, (d[str(i.lower()+j)].colour + 'o'), markerfacecolor='white')
#plt.xlim(0,100)
#plt.ylim(0,1)
plt.ylabel('% black in threshold')
plt.xlabel('PCBM concentration (%)')
#%%
'''
plot ratios of threshold images values from earlier defined python analysis function
'''
f = plt.figure(figsize=(4,4))
for i in experiments:
    x = []
    for j in ratios:
        x.append(d[str(i.lower()+j)].threshV)
    plt.plot(ratios, x, 'k-')
plt.xlim(0,100)
plt.ylim(0,1)
plt.ylabel('Black:white ratio')
plt.xlabel('PCBM concentration (%)')

#%%
'''
plotting hist
'''
#f, ax = plt.subplots(len(experiments), len(ratios))
i=1
f = plt.figure(figsize=(10,8))
for j in d.keys():
    ax = plt.subplot(len(experiments), len(ratios), i)
    d[j].hist.plot.bar(color='k')
    ax.set_xlim(0,200)
    ax.set_ylim(0,1E6)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticklabels([0,100])
    if (i-1)%5 != 0:
        ax.yaxis.set_ticklabels([])
    if i <15:
        ax.xaxis.set_ticklabels([])
    i+=1
f.subplots_adjust(wspace=0, hspace=0)

#%% plotting images
i = 1
f = plt.figure(figsize=(10,8))
for j in d.keys():
    ax = plt.subplot(len(experiments), len(ratios), i)
    if not d[j].imbw == None:
        ax.imshow(d[j].imbw, 'Greys')
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticklabels([])
    if (i-1)%5 != 0:
        ax.yaxis.set_ticklabels([])
    if i <15:
        ax.xaxis.set_ticklabels([])
    i+=1
f.subplots_adjust(wspace=0.05, hspace=0.05)

        
#%% plotting fouriers
saveDir = path + 'plots/fft/'
try:
    os.makedirs(saveDir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
        
for i in d.keys():
    F2 = d[i].fourier
    im = d[i].im
    r = d[i].fourierProfile
    middle = [F2.shape[0] / 2, F2.shape[1] / 2]
    a = 1.8 * len(F2)**2 
    b = 60 / len(F2)
    f, (ax1, ax3) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]}, figsize = (3,5))
    f.subplots_adjust(hspace=0)
    ax1.imshow(im)
    ax1.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
    ax1.set_title(d[i].name)
    ax1.axis('off')
    ax2 = f.add_axes([0.6, 0.61, 0.25, 0.25])
    ax2.imshow(np.abs(F2), 'gray', vmax= a)
    ax2.set_ylim((1-b)*middle[0],(1+b)*middle[0])
    ax2.set_xlim((1-b)*middle[1],(1+b)*middle[1])
    ax2.tick_params(axis = 'both', bottom = 'off', left = 'off', labelbottom='off', labelleft='off')
    ax2.axis('off')
    ax3.plot(r, 'ko-', markerfacecolor='white')
    ax3.set_ylim(-0.1*a, a)
    ax3.set_xlim(-0.02*(b * 2 * middle[0]), b * 2 * middle[0])
    ax3.set_yticks([])
    f.savefig(saveDir + i, bbox_inches='tight',dpi=300)       
        

#%%
data = pd.read_csv('H:/Research/OPV/Neutron Scattering/AFM/120degCheating/D80/cropped/D80_0min2.txt', encoding = 'ANSI')
#data.columns = [name]
#data = a['Height(Nanometer)\t'][45:].str.strip(to_strip='\t')
#data = pd.to_numeric(a)
split = np.array([data['Height(Nanometer)\t'].values.tolist()[i:i + 256] for i in range(0, len(data), 256)])

#%% create image
f, ax = plt.subplots()
ax.imshow(split, 'afmhot', vmin = -10, vmax = 10)

crossSectionPoint = 128# int(input('pixel of crossSeciton: '))

ax.plot([0,256],[crossSectionPoint,crossSectionPoint], 'w')

plt.plot(split[crossSectionPoint], 'k')

#%% create histogram
hist = np.histogram(data, bins = np.linspace(-10, 10, 100), density = True)[0]
plt.plot(hist)

#%% threshold
data2 = np.array(split)
data2[data2 <= 0] = 0
data2[data2 > 0] = 255
ratio = data2[data2>0].size/data2.size

#%% roughness
roughnessRMS = np.sqrt(np.mean(np.square(data.values)))         
roughnessRq = np.mean((np.abs(data.values)/len(data)))
