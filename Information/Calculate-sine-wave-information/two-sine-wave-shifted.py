import numpy as np
import matplotlib.pyplot as plt

def sine_wave(Steps=1000, Final=2*np.pi, Noise_level=0.1,Phase_Shift=np.pi/6):
    x = np.linspace(0, Final, Steps+1)
    y = np.sin(x+Phase_Shift) + Noise_level * np.random.randn(Steps+1)
    return x, y

def Frequency_wave(x,y):
    xs=x[2]-x[0]
    yf=[]
    xf=[]
    for i in range(len(y)-2):
        yf.append((y[i+2]-y[i])/xs)
        xf.append(x[i+1])
    return xf,yf

def Wave_Slop(y):
    xs=[i for i in range (len(y))]
    return xs,y

def Output_file(x,y,Name_file='Output',Format='txt'):
    with open(Name_file+'.'+Format, 'w') as file:
        for i in range(len(x)):
            file.write(f'{x[i]}\t{y[i]}\n')
    pass

def plot(x1,y1,x2,y2,address_name):
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.axhline(0, color = 'r', linestyle = '-.')
    #plt.show()
    plt.subplots_adjust(top = 0.98, bottom=0.08,left=0.06,right=0.98, hspace=0.2, wspace=0.44)
    plt.gcf().set_size_inches(12, 6.5)# don't change it
    plt.savefig(address_name+'.png',dpi=300)
    #plt.close()
    pass

def main():
    address_name='./input_data/two_sine_wave_with_constant_noise'
    x1,y1=sine_wave(10000,20*np.pi,0.1,0)
    #x2,y2=sine_wave(10000,20*np.pi,0,0)
    constant=0.7
    yc=[constant for i in range (len(x1))]

    for i in range (len(yc)):
        yc[i] = yc[i] + 0.1 * np.random.random() #np.random.rand(1)
    print(yc)


        #Output_file(x1,y1,'./Datas/sin1')
        #Output_file(x2,y2,'./Datas/sin2')
    Output_file(y1,yc,address_name)
        #x1f,y1f=Frequency_wave(x1,y1)
        #Output_file(x1f,y1f,'sin1f')
        #x1s,y1s=Wave_Slop(x1)
        #Output_file(x1s,y1s,'sin1s')
        #x2s,y2s=Wave_Slop(x2+np.pi/6)
        #Output_file(x2s,y2s,'sin2s')
    plot(x1,y1,x1,yc,address_name)
    print("Done :)")
    pass

if __name__=="__main__":
    main()