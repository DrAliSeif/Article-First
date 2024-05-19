import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
import math
import cmath
from scipy.signal import hilbert
from matplotlib.pyplot import figure


in1=0.01
out1=0.045

in2=0.08
out2=0.18

in3=0.5
out3=1.0

font1 = {'family': 'serif', 'color': 'blue', 'size': 190}
font2 = {'family': 'serif', 'color': 'darkred', 'size': 12}
fonts = 'Times New Roman'
dataset = np.genfromtxt(fname='right_sync.txt',skip_header=1)
xsignal = dataset[:,0]
ysignal = dataset[:,1]
signal_shift=ysignal

#t = np.linspace(0, 10, 40000)
time = np.linspace(0, 400,40000)

# Compute the Fourier Transform
freq = np.fft.fftfreq(len(signal_shift), time[1]-time[0])
ft = np.fft.fft(signal_shift)

#ft[0]=ft[1]
freq=freq#*1.89


fig = plt.figure()
#((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((signal))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
ax_1 = plt.subplot(5, 4, (1,4))
plt.plot(xsignal, signal_shift, color="#2F2CBF",linewidth = '2.8')
#plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, pad=15)
#plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=15)
#plt.xticks(np.arange(100, 385, 25),font=fonts,fontsize=31,color= '#262626')
#plt.yticks(font=fonts,fontsize=31,color= '#262626')
#plt.patch.set_visible(False)
plt.xlim(0,1)
plt.xlim(0,400)
plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')
plt.ylabel('Syncroney (r)',  fontsize=38, labelpad=37)
plt.xlabel('Time (s)',  fontsize=38)
#plt.xticks([])
#plt.yticks([])
#plt.box(False)
#((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((fft))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))




ax_2 = plt.subplot(5, 4, (5,8))
fft_vals = np.fft.fft(ysignal)
sampling_freq = 1 / (xsignal[1] - xsignal[0])  # Assuming evenly spaced samples
freqs = np.fft.fftfreq(len(ysignal), 1 / sampling_freq)
#plt.plot(freqs,  savgol_filter(np.abs(fft_vals), 80, 2), label='Fourier Transform',color='r')
fft_vals[0]=0
plt.plot(freqs,  savgol_filter(np.abs(fft_vals), 15, 2), label='Filtered Fourier Transform',color='#000077',linewidth = '4.4')
plt.fill(freqs,  savgol_filter(np.abs(fft_vals), 15, 2),color='#3C3CC7',alpha=0.8)
plt.fill_betweenx([-100, 100000], in2, out2, color="#007600", alpha=0.5)
plt.fill_betweenx([-100, 100000], in1, out1, color="#262728", alpha=0.5)


plt.fill_betweenx([-100, 100000], in3, out3, color="#FF3232", alpha=0.5)


plt.axvline(x=in1, color="#000000", alpha=0.7)
plt.axvline(x=(in1+out1)/2, color="#000000", linestyle='--', alpha=0.7)
plt.axvline(x=out1, color="#000000", alpha=0.7)


plt.axvline(x=in2, color="#003A00", alpha=0.7)
plt.axvline(x=(in2+out2)/2, color="#003A00", linestyle='--', alpha=0.7)
plt.axvline(x=out2, color="#003A00", alpha=0.7)

plt.axvline(x=in3, color="#800000", alpha=0.7)
plt.axvline(x=(in3+out3)/2, color="#800000", linestyle='--', alpha=0.7)
plt.axvline(x=out3, color="#800000", alpha=0.7)


'''plt.axvline(x=in1, color='k', linestyle='--', alpha=0.7)
plt.axvline(x=out1, color='k', linestyle='--', alpha=0.7)
plt.axvline(x=in2, color='k', linestyle='--', alpha=0.7)
plt.axvline(x=0.15, color='k', linestyle='--', alpha=0.7)'''
'''plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, pad=15)
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=15)
plt.xticks(font=fonts,fontsize=31,color= '#262626')
plt.yticks(font=fonts,fontsize=31,color= '#262626')
plt.fill_betweenx([-100, 100000], in1, out1, color="#262728", alpha=0.5)
plt.fill_betweenx([-100, 100000], out1, in2, color="b", alpha=0.5)
plt.fill_betweenx([-100, 100000], in2, out2, color="#FFFF00", alpha=0.5)
plt.fill_betweenx([-100, 100000], in3, out3, color="g", alpha=0.5)'''


'''plt.axvline(x=out1, color='k', linestyle='dashed')
plt.axvline(x=0.15, color='k', linestyle='dashed')
plt.axvline(x=in2, color='k', linestyle='dashed')'''


plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')

plt.ylabel('Syncroney (r)',  fontsize=38, labelpad=37)
plt.xlabel('Time (s)',  fontsize=38)
#plt.box(False)
plt.ylim(0, 1700)
#plt.xlim(0, 4)
plt.xlim(0, 4)
#plt.xticks(np.arange(0, 1, 0.05),font=fonts,fontsize=31,color= '#262626')
#plt.yticks([])
#plt.xticks([])


'''ax2 = fig.add_axes([0.5, 0.68, 0.4, 0.1])#{x,y,l,h}
ax2.plot(time, savgol_filter(np.abs(fft_vals), 15, 2), color='#000077', label='cos(angle(H(real(Ifft(fft(xR),in2-0.15))))-angle(H(real(Ifft(fft(xL,xR),in2-0.15)))))',alpha=1)
ax2.fill_between(time, savgol_filter(np.abs(fft_vals), 15, 2), -1,color='#3C3CC7',alpha=0.6)
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1900)

ax2.fill_betweenx([-100, 100000], in1, out1, color="#262728", alpha=0.5)
#ax2.fill_betweenx([-100, 100000], out1, in2, color="b", alpha=0.5)
ax2.fill_betweenx([-100, 100000], in2, out2, color="#007600", alpha=0.5)'''

#ax2.axis('off')



#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))





ax_3 = plt.subplot(5, 4, (9,12))
#______________________________________________________filter     
'''fill_startr = 2.3
fill_endr = out3
fill_midr = 3
freq_min_r = fill_startr  # Hz
freq_max_r = fill_endr  # Hz
filter_r = np.logical_or(freqs < freq_min_r, freqs > freq_max_r)
fft_vals_filtered_r = fft_vals.copy()
fft_vals_filtered_r[filter_r] = 0
signal_filtered_r = np.real(np.fft.ifft(fft_vals_filtered_r))
plt.plot(xsignal, signal_filtered_r, 'r', label='Signal(Red)')
#ax[5].plot(xsignal, ysignal-(np.average(ysignal)), '#262626', label='Signal(x)',alpha=0.3)
plt.legend(loc="upper right",prop={'family': 'Times New Roman',  "size": 20 })
plt.xlim(0, 400)
#ax[5] = ax[5].twinx()'''

#______________________________________________________filter     
freq_min_g = in1  # Hz
freq_max_g = out1  # Hz
filter_g = np.logical_or(freqs < freq_min_g, freqs > freq_max_g)
fft_vals_filtered_g = fft_vals.copy()
fft_vals_filtered_g[filter_g] = 0
signal_filtered_g = np.real(np.fft.ifft(fft_vals_filtered_g))
plt.plot(xsignal, signal_filtered_g, color='#262728', label='Signal(Green)',linewidth = '4.4')




from scipy.signal import hilbert, chirp
phase_A = np.angle(hilbert(signal_filtered_g))
amps_B = np.abs(hilbert(signal_filtered_g))#SATR
ffp = phase_A/20
AfA = amps_B



plt.plot(xsignal, ffp, color='#00FFFF', label='Signal(Green)', linestyle='dashed',linewidth = '4.4')
plt.plot(xsignal, AfA, color='#873200', label='Signal(Green)', linestyle='dashed',linewidth = '4.4')

#ax[5].plot(xsignal, ysignal-(np.average(ysignal)), '#262626', label='Signal(x)',alpha=0.3)
#plt.legend(loc="upper right",prop={'family': 'Times New Roman',  "size": 20 })
plt.xlim(0, 400)

plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')

plt.ylabel('Syncroney (r)',  fontsize=38, labelpad=37)
plt.xlabel('Time (s)',  fontsize=38)
#plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, pad=15)
#plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=15)
#plt.xticks(font=fonts,fontsize=31,color= '#262626')
#plt.yticks(font=fonts,fontsize=31,color= '#262626')
#plt.box(False)
#plt.yticks([])
#plt.xticks([])
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


#ax_4 = plt.subplot(5, 4, (13,16))
'''freq_min_g = in1  # Hz
freq_max_g = out1  # Hz
filter_g = np.logical_or(freqs < freq_min_g, freqs > freq_max_g)
fft_vals_filtered_g = fft_vals.copy()
fft_vals_filtered_g[filter_g] = 0
signal_filtered_g = np.real(np.fft.ifft(fft_vals_filtered_g))
plt.plot(xsignal, signal_filtered_g, 'k', label='Signal(Green)')
plt.plot(xsignal, AfA, 'r', label='Signal(Green)')

#ax[5].plot(xsignal, ysignal-(np.average(ysignal)), '#262626', label='Signal(x)',alpha=0.3)
#plt.legend(loc="upper right",prop={'family': 'Times New Roman',  "size": 20 })
plt.xlim(0, 400)

plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')'''

#plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, pad=15)
#plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=15)
#plt.xticks(font=fonts,fontsize=31,color= '#262626')
#plt.yticks(font=fonts,fontsize=31,color= '#262626')
#plt.box(False)
#plt.yticks([])
#plt.xticks([])

#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
#((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((top polar))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

ax_5 = plt.subplot(5, 4, 16, projection='polar')
#______________________________________________________filter      _A
freq_min_A = in1  # Hz
freq_max_A = out1  # Hz
filter_A = np.logical_or(freq < freq_min_A, freq > freq_max_A)
fft_vals_filtered_A = ft.copy()
fft_vals_filtered_A[filter_A] = 0
signal_filtered_A = np.real(np.fft.ifft(fft_vals_filtered_A))
#______________________________________________________filter      _B
freq_min_B = in2  # Hz0.06, 0.15
freq_max_B = out2 # Hz
filter_B = np.logical_or(freq < freq_min_B, freq > freq_max_B)
fft_vals_filtered_B = ft.copy()
fft_vals_filtered_B[filter_B] = 0
signal_filtered_B = np.real(np.fft.ifft(fft_vals_filtered_B))


phase_A = np.angle(hilbert(signal_filtered_A))

amps_B = np.abs(hilbert(signal_filtered_B))#SATR
phase_delta = np.angle(hilbert(signal_filtered_A))
amps_gamma = np.abs(hilbert(signal_filtered_B))#SATR
phase_amps_gamma= np.angle(hilbert(amps_gamma))
defrent_each_angles=phase_delta-phase_amps_gamma



#calculate Average radius
def order_parameter(N, phi):
    rc = 0.0
    rs = 0.0
    for j in range(N):
        rc += math.cos(phi[j])
        rs += math.sin(phi[j])
    return math.sqrt(rc**2 + rs**2) / (1.0 * N)
average_radius2=order_parameter(len(defrent_each_angles),defrent_each_angles)



#calculate Average angels
resolt_rad=0
average_angels=0
for i in range(0,len(defrent_each_angles)):
    resolt_rad=resolt_rad+cmath.exp(complex(0, defrent_each_angles[i]))
resolt_rad=resolt_rad/len(defrent_each_angles)
if resolt_rad.real<0:
    resolt_rad=cmath.atan(resolt_rad.imag/resolt_rad.real)#resolt_degree=resolt_rad*180/math.pi
    average_angels=resolt_rad.real+math.pi
else:
    resolt_rad=cmath.atan(resolt_rad.imag/resolt_rad.real)#resolt_degree=resolt_rad*180/math.pi
    average_angels=resolt_rad.real


print("PLV out2=")
print(average_radius2)

plt.plot([defrent_each_angles,defrent_each_angles], [0,1],alpha=0.005, color='k')
plt.plot([average_angels,average_angels], [0,np.abs(average_radius2)],alpha=1, color='r', linewidth = '4')    
plt.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, labeltop=False,pad=26,  labelsize=32, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=32, labelcolor='#262626')
plt.yticks([])
#plt.xticks([])
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
ax_6 = plt.subplot(5, 4, 15)
from scipy import stats
#______________________________________________________filter      _A
freq_min_A = in1  # Hz
freq_max_A = out1  # Hz
filter_A = np.logical_or(freq < freq_min_A, freq > freq_max_A)
fft_vals_filtered_A = ft.copy()
fft_vals_filtered_A[filter_A] = 0
signal_filtered_A = np.real(np.fft.ifft(fft_vals_filtered_A))
#______________________________________________________filter      _B
freq_min_B = in2  # Hz
freq_max_B = out2 # Hz
filter_B = np.logical_or(freq < freq_min_B, freq > freq_max_B)
fft_vals_filtered_B = ft.copy()
fft_vals_filtered_B[filter_B] = 0
signal_filtered_B = np.real(np.fft.ifft(fft_vals_filtered_B))


phase_A = np.angle(hilbert(signal_filtered_A))

amps_B = np.abs(hilbert(signal_filtered_B))#SATR
# Bin the phases and calculate the mean amplitude for each bin
phase_bins = np.linspace(-np.pi, np.pi, 18) # example number of bins
kAfAlffp, _, _ = stats.binned_statistic(ffp, AfA, bins=phase_bins, statistic='mean')
n, bins = np.histogram(phase_A, bins=17, range=(-np.pi,np.pi), weights=amps_B, density=False)
P = kAfAlffp / np.sum(kAfAlffp)
P2 = n / np.sum(n)
# Calculate the MI
MI_2=0
for i in range (0,len(P2)):
    if P2[i]!=0:
        MI_2 = MI_2+(P2[i] * np.log10(P2[i]))
    else:
        MI_2 = MI_2
MI_2=1+(MI_2 / np.log10(len(P2)))
print("MI out2=")
print(MI_2)
plt.bar(phase_bins[:-1], P2,edgecolor='k', width=.366, color = u'#7C7C7C')


plt.xlim(-np.pi-0.19, np.pi-0.19)
plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')
plt.ylabel('Syncroney (r)',  fontsize=38)
plt.xlabel('Time (s)',  fontsize=38)
#plt.box(False)


#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

ax_7 = plt.subplot(5, 4,  (13,14))

freq_min_g = in1  # Hz
freq_max_g = out1  # Hz
filter_g = np.logical_or(freqs < freq_min_g, freqs > freq_max_g)
fft_vals_filtered_g = fft_vals.copy()
fft_vals_filtered_g[filter_g] = 0
signal_filtered_g = np.real(np.fft.ifft(fft_vals_filtered_g))
plt.plot(xsignal, signal_filtered_g, color='#262728', label='Signal(Green)',linewidth = '4')

#plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, pad=15)
#plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=15)
#plt.xticks(font=fonts,fontsize=31,color= '#262626')
#plt.yticks(font=fonts,fontsize=31,color= '#262626')
#plt.box(False)
#plt.yticks([])
#plt.xticks([])
freq_min_g = in2  # Hz
freq_max_g = out2 # Hz
filter_g = np.logical_or(freqs < freq_min_g, freqs > freq_max_g)
fft_vals_filtered_g = fft_vals.copy()
fft_vals_filtered_g[filter_g] = 0
signal_filtered_g = np.real(np.fft.ifft(fft_vals_filtered_g))
plt.plot(xsignal, signal_filtered_g, color= '#007600', label='Signal(Green)',linewidth = '4')

phase_A = np.angle(hilbert(signal_filtered_g))
amps_B = np.abs(hilbert(signal_filtered_g))#SATR
ffp = phase_A/20
AfA = amps_B



#plt.plot(xsignal, ffp, color='#00FFFF', label='Signal(Green)', linestyle='dashed',linewidth = '4.4')
plt.plot(xsignal, AfA, color='#873200', label='Signal(Green)', linestyle='dashed',linewidth = '4.4')


plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')
plt.xlim(0, 400)
#plt.ylim(0, 1)
plt.ylabel('Syncroney (r)',  fontsize=38, labelpad=37)
plt.xlabel('Time (s)',  fontsize=38)
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
#((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((2))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

ax_8 = plt.subplot(5, 4,  20, projection='polar')

#______________________________________________________filter      _A
freq_min_A = in1  # Hz
freq_max_A = out1  # Hz
filter_A = np.logical_or(freq < freq_min_A, freq > freq_max_A)
fft_vals_filtered_A = ft.copy()
fft_vals_filtered_A[filter_A] = 0
signal_filtered_A = np.real(np.fft.ifft(fft_vals_filtered_A))
#______________________________________________________filter      _B
freq_min_B = in3  # Hz
freq_max_B = out3 # Hz
filter_B = np.logical_or(freq < freq_min_B, freq > freq_max_B)
fft_vals_filtered_B = ft.copy()
fft_vals_filtered_B[filter_B] = 0
signal_filtered_B = np.real(np.fft.ifft(fft_vals_filtered_B))


phase_A = np.angle(hilbert(signal_filtered_A))
amps_B = np.abs(hilbert(signal_filtered_B))#SATR

phase_delta = np.angle(hilbert(signal_filtered_A))
amps_gamma = np.abs(hilbert(signal_filtered_B))#SATR

phase_amps_gamma= np.angle(hilbert(amps_gamma))
defrent_each_angles=phase_delta-phase_amps_gamma



#calculate Average radius
def order_parameter(N, phi):
    rc = 0.0
    rs = 0.0
    for j in range(N):
        rc += math.cos(phi[j])
        rs += math.sin(phi[j])
    return math.sqrt(rc**2 + rs**2) / (1.0 * N)
average_radius=order_parameter(len(defrent_each_angles),defrent_each_angles)

print("PLV out3=")
print(average_radius)

#calculate Average angels
resolt_rad=0
average_angels=0
for i in range(0,len(defrent_each_angles)):
    resolt_rad=resolt_rad+cmath.exp(complex(0, defrent_each_angles[i]))
resolt_rad=resolt_rad/len(defrent_each_angles)
if resolt_rad.real<0:
    resolt_rad=cmath.atan(resolt_rad.imag/resolt_rad.real)#resolt_degree=resolt_rad*180/math.pi
    average_angels=resolt_rad.real+math.pi
else:
    resolt_rad=cmath.atan(resolt_rad.imag/resolt_rad.real)#resolt_degree=resolt_rad*180/math.pi
    average_angels=resolt_rad.real
    



plt.plot([defrent_each_angles,defrent_each_angles], [0,1],alpha=0.005, color='k')
plt.plot([average_angels,average_angels], [0,np.abs(average_radius)],alpha=1, color='r', linewidth = '4')    
plt.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, labeltop=False,pad=26,  labelsize=32, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=32, labelcolor='#262626')

plt.yticks([])
#plt.xticks([])
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

ax_9 = plt.subplot(5, 4,  19)
plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, pad=15)
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=15)

#______________________________________________________filter      _A
freq_min_A = in1  # Hz
freq_max_A = out1  # Hz
filter_A = np.logical_or(freq < freq_min_A, freq > freq_max_A)
fft_vals_filtered_A = ft.copy()
fft_vals_filtered_A[filter_A] = 0
signal_filtered_A = np.real(np.fft.ifft(fft_vals_filtered_A))
#______________________________________________________filter      _B
freq_min_B = in3  # Hz
freq_max_B = out3 # Hz
filter_B = np.logical_or(freq < freq_min_B, freq > freq_max_B)
fft_vals_filtered_B = ft.copy()
fft_vals_filtered_B[filter_B] = 0
signal_filtered_B = np.real(np.fft.ifft(fft_vals_filtered_B))


phase_A = np.angle(hilbert(signal_filtered_A))

amps_B = np.abs(hilbert(signal_filtered_B))#SATR
# Bin the phases and calculate the mean amplitude for each bin
phase_bins = np.linspace(-np.pi, np.pi, 18) # example number of bins
kAfAlffp, _, _ = stats.binned_statistic(ffp, AfA, bins=phase_bins, statistic='mean')
n, bins = np.histogram(phase_A, bins=17, range=(-np.pi,np.pi), weights=amps_B, density=False)
P = kAfAlffp / np.sum(kAfAlffp)
P2 = n / np.sum(n)


# Calculate the MI
MI_1=0
for i in range (0,len(P2)):
    if P2[i]!=0:
        MI_1 = MI_1+(P2[i] * np.log10(P2[i]))
    else:
        MI_1 = MI_1

MI_1=1+(MI_1 / np.log10(len(P2)))


print("MI out3=")
print(MI_1)
plt.bar(phase_bins[:-1], P2,edgecolor='k', width=.366, color = '#7C7C7C')
#plt.box(False)
plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')
plt.xlim(-np.pi-0.19, np.pi-0.19)
#plt.ylim(0, 1)
plt.ylabel('Syncroney (r)',  fontsize=38)
plt.xlabel('Time (s)',  fontsize=38)
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

ax_10 = plt.subplot(5, 4,  (17,18))

freq_min_a = in1  # Hz
freq_max_a = out1  # Hz
filter_a = np.logical_or(freqs < freq_min_a, freqs > freq_max_a)
fft_vals_filtered_a = fft_vals.copy()
fft_vals_filtered_a[filter_a] = 0
signal_filtered_a = np.real(np.fft.ifft(fft_vals_filtered_a))
plt.plot(xsignal, signal_filtered_a, color='#262728',linewidth = '4')


freq_min_b = in3  # Hz
freq_max_b = out3 # Hz
filter_b = np.logical_or(freqs < freq_min_b, freqs > freq_max_b)
fft_vals_filtered_b = fft_vals.copy()
fft_vals_filtered_b[filter_b] = 0
signal_filtered_b = np.real(np.fft.ifft(fft_vals_filtered_b))
plt.plot(xsignal, signal_filtered_b, color='#FF3232',linewidth = '4')

'''phase_A = np.angle(hilbert(signal_filtered_g))
amps_B = np.abs(hilbert(signal_filtered_g))#SATR
ffp = phase_A/20
AfA = amps_B
#plt.plot(xsignal, ffp, color='#00FFFF', label='Signal(Green)', linestyle='dashed',linewidth = '4.4')
#plt.plot(xsignal, AfA, color='#873200', label='Signal(Green)', linestyle='dashed',linewidth = '4.4')'''

plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,  labelsize=37, labelcolor='#262626')
plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False,  labelsize=37, labelcolor='#262626')
plt.xlim(0, 400)
#plt.ylim(0, 1)
plt.ylabel('Syncroney (r)',  fontsize=38, labelpad=37)
plt.xlabel('Time (s)',  fontsize=38)
#(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
'''ax2 = fig.add_axes([0.5, 0.68, 0.4, 0.1])#{x,y,l,h}
#ax2.plot(time,signal_shift-0.5 , color='b',alpha=1)
ax2.plot(time,0.5+(signal_filtered_b+signal_filtered_g+signal_filtered_a)*2 , color='r',alpha=1)
ax2.set_xlim(0, 400)'''

plt.subplots_adjust(right=0.91,left= 0.07,top = 0.98, bottom=0.08, hspace=0.38, wspace=0.59)
plt.gcf().set_size_inches(38, 28)

name=f"{in1:.3f}-{out1:.3f}-{in2:.3f}-{out2:.3f}-{in3:.3f}-{out3:.3f}-PVL({average_radius2:.3f},{average_radius:.3f})-MI({MI_2:.3f},{MI_1:.3f})"#-str(round(out3,4))+"=PVL"+str(average_radius)+"-MI"+str(MI_1)
plt.savefig('./'+name+'.png', dpi=100, bbox_inches='tight', pad_inches=1, bbox_extra_artists=[])
