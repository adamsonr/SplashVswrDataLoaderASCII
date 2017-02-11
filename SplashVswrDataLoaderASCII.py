import numpy as np

from pylab import *



class pyVswrSourceData:
    def __init__(self, asciiData):
        
        #Retain a copy of the c-struct data object        
        self.asciiData = asciiData
        
   
        #Count the total number of recods present    
        recordCount = asciiData.shape[0]
        self.recordCount = recordCount
        
        #Allocate numpy arrays to receive data
        self.nPhases = 6
        self.allocateNumpy()
        
        #Copy values from C-Struct to Numpy Arrays
        self.convertToNumpy()
        
        
    def allocateNumpy(self):
        
        recordCount = self.recordCount
        
        self.tcpId = np.zeros(recordCount,dtype=complex128)
        self.recordId = np.zeros(recordCount,dtype=complex128)

        #DC Offset Values on raw measurements
        #Alias provided for compatibility with existing processing routines
        self.Ioff_i = self.offsetCurrentI = np.zeros(recordCount,dtype=complex128)
        self.Ioff_q = self.offsetCurrentQ = np.zeros(recordCount,dtype=complex128)
        self.Voff_i = self.offsetVoltageI = np.zeros(recordCount,dtype=complex128)
        self.Voff_q = self.offsetVoltageQ = np.zeros(recordCount,dtype=complex128)
        
        #Raw IQ Values for voltage and current
        #Phase values are combine to reconstruct harmonic content
        N_RECORDS = self.recordCount
        N_PHASES = self.nPhases
        self.I_i = np.zeros([N_RECORDS,N_PHASES],dtype=complex128)
        self.I_q = np.zeros([N_RECORDS,N_PHASES],dtype=complex128)
        self.V_i = np.zeros([N_RECORDS,N_PHASES],dtype=complex128)
        self.V_q = np.zeros([N_RECORDS,N_PHASES],dtype=complex128)        
        
        #Alias for compatibility with existing processing routines        
        self.raw00CurrentI = self.I_i[:,0]
        self.raw00CurrentQ = self.I_q[:,0]
        self.raw00VoltageI = self.V_i[:,0]
        self.raw00VoltageQ = self.V_q[:,0]
        
        self.raw15CurrentI = self.I_i[:,1]
        self.raw15CurrentQ = self.I_q[:,1]
        self.raw15VoltageI = self.V_i[:,1]
        self.raw15VoltageQ = self.V_q[:,1]
        
        self.raw30CurrentI = self.I_i[:,2]
        self.raw30CurrentQ = self.I_q[:,2]
        self.raw30VoltageI = self.V_i[:,2]
        self.raw30VoltageQ = self.V_q[:,2]
        
        self.raw45CurrentI = self.I_i[:,3]
        self.raw45CurrentQ = self.I_q[:,3]
        self.raw45VoltageI = self.V_i[:,3]
        self.raw45VoltageQ = self.V_q[:,3]
        
        self.raw60CurrentI = self.I_i[:,4]
        self.raw60CurrentQ = self.I_q[:,4]
        self.raw60VoltageI = self.V_i[:,4]
        self.raw60VoltageQ = self.V_q[:,4]
        
        self.raw75CurrentI = self.I_i[:,5]
        self.raw75CurrentQ = self.I_q[:,5]
        self.raw75VoltageI = self.V_i[:,5]
        self.raw75VoltageQ = self.V_q[:,5]

        

    def convertToNumpy(self):
        
        
        
        self.tcpId = self.asciiData[:,1]
        self.recordId = self.asciiData[:,2]
        self.offsetCurrentI[:] = self.asciiData[:,3]
        self.offsetCurrentQ[:] = self.asciiData[:,4]
        self.offsetVoltageI[:] = self.asciiData[:,5]
        self.offsetVoltageQ[:] = self.asciiData[:,6]
        
        self.raw00CurrentI[:] = self.asciiData[:,7]
        self.raw00CurrentQ[:] = self.asciiData[:,8]
        self.raw00VoltageI[:] = self.asciiData[:,9]
        self.raw00VoltageQ[:] = self.asciiData[:,10]
        
        self.raw15CurrentI[:] = self.asciiData[:,11]
        self.raw15CurrentQ[:] = self.asciiData[:,12]
        self.raw15VoltageI[:] = self.asciiData[:,13]
        self.raw15VoltageQ[:] = self.asciiData[:,14]
        
        self.raw30CurrentI[:] = self.asciiData[:,15]
        self.raw30CurrentQ[:] = self.asciiData[:,16]
        self.raw30VoltageI[:] = self.asciiData[:,17]
        self.raw30VoltageQ[:] = self.asciiData[:,18]
        
        self.raw45CurrentI[:] = self.asciiData[:,19]
        self.raw45CurrentQ[:] = self.asciiData[:,20]
        self.raw45VoltageI[:] = self.asciiData[:,21]
        self.raw45VoltageQ[:] = self.asciiData[:,22]
        
        self.raw60CurrentI[:] = self.asciiData[:,23]
        self.raw60CurrentQ[:] = self.asciiData[:,24]
        self.raw60VoltageI[:] = self.asciiData[:,25]
        self.raw60VoltageQ[:] = self.asciiData[:,26]
        
        self.raw75CurrentI[:] = self.asciiData[:,27]
        self.raw75CurrentQ[:] = self.asciiData[:,28]
        self.raw75VoltageI[:] = self.asciiData[:,29]
        self.raw75VoltageQ[:] = self.asciiData[:,30]


class pyVswrData:
    def __init__(self, asciiData,f0,reordering=(0,1,2,3,4,5,6,7,8,9,10,11)):
        
        #Assign constants
        self.nPhases = 6
        self.phases = np.linspace(0,90,self.nPhases,endpoint=False)
        
        #Unpack generic array data into named fields        
        self.sourceData = pyVswrSourceData(asciiData)
        
        #Bring a few of the frequently used itesm to top level
        self.recordCount = self.sourceData.recordCount
        
        self.tcpId = self.sourceData.tcpId
        self.recordId = self.sourceData.recordId
        self.f0=f0
        self.reordering=reordering
        
    def processVI(self,f0,
                  adcScale=5.12,
                  adcBits=16, 
                  RSense=0.1,
                  RF=3320.0, 
                  G2=0.001, 
                  G3=0.00217, 
                  Vff=1.,
                  Lcorr=0.,
                  Ccorr=0.,
                  tdelay=0.,
                  apply_correction=False):
        
        N_RECORDS = self.recordCount
        N_PHASES = self.nPhases
        
        self.vmeas=zeros([self.recordCount,2*N_PHASES])
        self.imeas=zeros([self.recordCount,2*N_PHASES])
        self.vreordered=zeros([self.recordCount,2*N_PHASES])
        self.ireordered=zeros([self.recordCount,2*N_PHASES])
        
        V = np.zeros([N_RECORDS,N_PHASES],dtype=np.complex128)
        I = np.zeros([N_RECORDS,N_PHASES],dtype=np.complex128)
        
        Voff = np.zeros(N_RECORDS,dtype=np.complex128)
        Ioff = np.zeros(N_RECORDS,dtype=np.complex128)
        
        #Translate I/Q to complex values        
        V = self.sourceData.V_i + 1j*self.sourceData.V_q
        I = self.sourceData.I_i + 1j*self.sourceData.I_q

        Voff = self.sourceData.Voff_i + 1j*self.sourceData.Voff_q
        Ioff = self.sourceData.Ioff_i + 1j*self.sourceData.Ioff_q
        
        #Subtract offsets
        for i in range(N_PHASES):
            V[:,i] -= Voff
            I[:,i] -= Ioff
        
        #Scale to get ADC output
        V = V*adcScale/(2.**adcBits)
        I = I*adcScale/(2.**adcBits)
        
        #Scale to get voltage and current
        self.RSense=RSense
        RSense=RSense   #Current sense resistance (Ohms)
        RF = RF         #transimpedance gain of TZA (Ohms)
        #G1 = RSense+1j*2*pi*1.7e-9*f0     #Gain of current sense resistor (Ohms)
        G1=RSense
        G2 = G2         #Gain of voltage divider for voltage measurement V/V
        G3 = G3         #S
        
        Vff = Vff       #Gain correction from comparison with scope measurements
        
        if apply_correction==True:
            V = V/(RF*G3*G2)*Vff
            I = I/(RF*G3*(G1+1j*f0*2*pi*Lcorr))
            I=I-V*(1j*2*pi*Ccorr*f0)
            I=I*exp(1j*2*pi*f0*tdelay)
        else:
            V = V/(RF*G3*G2)*Vff
            I = I/(RF*G3*G1)
        #Assign results as class members
        self.V = V
        self.I = I
        
    def processHarmonics(self, 
                         f0,
                         show_iq=False,
                         record_to_show=20,
                         which_meas=0):
        
        N_RECORDS = self.recordCount
        N_PHASES = self.nPhases
        
        V = self.V
        I = self.I
        
        #Harmonic Analysis
        #Stores the fundamental component of voltage and current
        V1 = np.zeros(N_RECORDS,dtype=np.complex128)
        I1 = np.zeros(N_RECORDS,dtype=np.complex128)
        #absz  = np.zeros(N_RECORDS,dtype=np.complex128)

        #Fourier components
        Ais,Avs = np.zeros([N_RECORDS,6],dtype=np.complex128),np.zeros([N_RECORDS,6],dtype=np.complex128)
        
        for N in range(N_RECORDS):
            Vf = V[N,:]
            If = I[N,:]
        
            #Generate the inversion matrix for 15 degree (pi/12 rad) steps
            F = np.zeros([12,12],dtype=np.float64)
               
            for ik,k in enumerate(np.arange(1,13,2)): #indexes the frequency
                for n in range(12): #Indexes the phase
                    #Re{A}
                    F[ik,n]=np.cos(np.pi/12*k*n)
                    #Im{A}
                    F[ik+6,n]=-np.sin(np.pi/12*k*n)
                F=np.matrix(F)
        
            #In Daxsonics lab collection scheme where phases shifts are applied
            #to the PA signal.  This code will have to be changed slightly when
            #phase shifts are applied to the LO signal.  Once I have some reference
            #data from BM I can figure out what the changes are
            m=np.matrix(np.concatenate((np.real(V[N,:]),-np.imag(V[N,:])))).T
            
            #Fudge to fix bad phase measurement ordering
            q=array(m)
            m=matrix([q[i] for i in self.reordering])
            self.vmeas[N,:]=q[:,which_meas]
            self.vreordered[N,:]=array([q[i,which_meas] for i in self.reordering])
            
        
            if show_iq==True:
                if N==record_to_show:
                    figure()
                    title("Voltage IQ")
                    plot(m)
                    grid('on')
                    xlabel("Measurement number")
                    ylabel("Voltage (V)")
                
            #Construct the Fourier amplitudes from the extracted measurements
            Av=(F*m)[0:6]+1j*(F*m)[6:12]  
            #Convert from matrix to array
            Av=np.array(Av.T)[0]/6*np.arange(1,13,2)*np.array([1,-1,1,-1,1,-1])
    
            #Calculate current Fourier amplitudes
            m=np.matrix(np.concatenate((np.real(If),-np.imag(If)))).T
            #Construct the Fourier amplitudes from the extracted measurements
            q=array(m)
            self.ireordered[N,:]=array([q[i,which_meas] for i in self.reordering])
            m=matrix([q[i] for i in self.reordering])
            if show_iq==True:
                if N==record_to_show:
                    figure()
                    title("Current IQ")
                    plot(m)
                    grid('on')
                    xlabel("Measurement number")
                    ylabel("Current (A)")
                    print "TCP id=",self.tcpId[N]
                    print "Record id=",self.recordId[N]
                    
            self.imeas[N,:]=q[:,which_meas]
            
            Ai=(F*m)[0:6]+1j*(F*m)[6:12]  
            #Convert from matrix to array
            Ai=np.array(Ai.T)[0]/6*np.arange(1,13,2)*np.array([1,-1,1,-1,1,-1])
            
            #Cshunt=360e-12
            #for q in range(6):
            #    Ai[q]=Ai[q]-Av[q]*1j*2*pi*Cshunt*(self.f0*(2*q+1))
        
            Ais[N,:]=Ai
            Avs[N,:]=Av
            
            V1[N]=Avs[N,0]
            I1[N]=Ais[N,0]
            
                
     
            
            #Assign results as class member
            self.V1 = V1
            self.I1 = I1
            self.Ais=Ais
            self.Avs=Avs
        show()
    
            
    def processOutputs(self, Z0=50.,apply_compensation=False):
        
        V1RMS = self.V1/sqrt(2)
        I1RMS = self.I1/sqrt(2)
       
        #load impedance @ fundamental        
        Zraw = V1RMS/I1RMS
        if apply_compensation:
            zcorr=self.CorrectZ(Zraw,f0)
        else:
            zcorr=Zraw
        I1RMS=V1RMS/zcorr
        #Load reflection coefficient
        RC=(V1RMS-Z0*I1RMS)/(V1RMS+Z0*I1RMS)
        
        #VSWR        
        VSWR=(1.+np.abs(RC))/(1.-np.abs(RC))
        
        #Forward and refelected power
        PF=np.abs(V1RMS+Z0*I1RMS)**2./(4.*Z0)
        PR=np.abs(V1RMS-Z0*I1RMS)**2./(4.*Z0)
        Pload=PF-PR
        #assign results as class members
        self.Zraw = Zraw
        self.Z=zcorr
        self.RC = RC
        self.VSWR = VSWR
        self.PF = PF
        self.PR = PR
        self.Pload=Pload
        
    def process_scope_traces(self,nrec,Npts=1000):
        '''
        Constructs arrays length Npts containing reconstructed current and voltage 
        waveforms
        '''
        self.scopephi=linspace(0,4*pi,Npts)
        self.scopev=sum(real(array([self.Avs[nrec,n]*exp(1j*(2*n+1)*self.scopephi)for n in range(len(self.Avs[nrec,:]-1))])),0)
        self.scopei=sum(real(array([self.Ais[nrec,n]*exp(1j*(2*n+1)*self.scopephi) for n in range(len(self.Ais[nrec,:]-1))])),0)
        
    def CorrectZ(self,zmeas,f0,cal_file="z_cal.csv"):
        '''
        Loads in calibration table and applies correction.  Cal table
        consists of rows of pre-calculated Z11, Z12, Z21 and Z22 values 
        '''
        #Load calibration table
        f0s,Z11sr,Z11si,Z12sr,Z12si,Z21sr,Z21si,Z22sr,Z22si=loadtxt(cal_file,delimiter=',')
        Z11s=Z11sr+1j*Z11si
        Z12s=Z12sr+1j*Z12si
        Z21s=Z21sr+1j*Z21si
        Z22s=Z22sr+1j*Z22si
        #Cubic spline interpolation
        from scipy.interpolate import interp1d
        Z11interp=interp1d(f0s,Z11s,kind='cubic')
        Z12interp=interp1d(f0s,Z12s,kind='cubic')
        Z22interp=interp1d(f0s,Z22s,kind='cubic')
        Z21interp=interp1d(f0s,Z21s,kind='cubic')
        #Calculate interpolated value of impedance parameters
        Z11=Z11interp(f0)
        Z12=Z12interp(f0)
        Z22=Z22interp(f0)
        Z21=Z21interp(f0)
        
        #Apply correction
        zcorr=Z12*Z21/(Z11-zmeas)-Z22
        return(zcorr)        

    
    
    
if __name__ == "__main__":
    #Close plots
    close('all')
    #Define constants
    use_tcps=(200,)    
    import os
    reorderings={}
    reorderings['10MHz']=[3,4,5,0,1,2,9,10,11,6,7,8]
    reorderings['11MHz']=[3,4,5,0,1,2,9,10,11,6,7,8]
    reorderings['12MHz']=[3,4,5,0,1,2,9,10,11,6,7,8]
    reorderings['2MHz']=[5,0,1,2,3,4,11,6,7,8,9,10]
    reorderings['3MHz']=[0,1,2,3,4,5,6,7,8,9,10,11]
    reorderings['4MHz']=[1,2,3,4,5,0,7,8,9,10,11,6]
    reorderings['5MHz']=[1,2,3,4,5,0,7,8,9,10,11,6]
    reorderings['6MHz']=[2,3,4,5,0,1,8,9,10,11,6,7]
    reorderings['7MHz']=[2,3,4,5,0,1,8,9,10,11,6,7]
    reorderings['8MHz']=[2,3,4,5,0,8,9,10,11,6,1,7]
    reorderings['9MHz']=[3,4,5,0,1,2,9,10,11,6,7,8]
    
    #reorderings['8MHz']=[0,1,2,3,4,5,6,7,8,9,10,11]
#    reorderings['9MHz']=[0,1,2,3,4,5,6,7,8,9,10,11]
    reorderings['11MHz']=[0,1,2,3,4,5,6,7,8,9,10,11]
    
    
    
    
    for (filedir,_,filenames) in os.walk('./'):
#        if filedir!=('./'): #ignore root-level files
        if filedir==('./10 Ohms'): #ignore root-level files
            for filename in filenames:
                print filename
                f0=float(filename.split('MHz')[0])*1e6
                fullname=filedir+'/'+filename
                #Load the file    
                monitorData = np.loadtxt(fullname,dtype=np.uint16,delimiter=",")
                freqstr=filename.split('_')[0]
                #Translate generic array data to named fields
                tData = pyVswrData(monitorData,f0,reordering=reorderings[freqstr])
    
                #Post process raw values
                tData.processVI(f0)
    
                #Calculate spectral components from phase data
                tData.processHarmonics(f0=f0,show_iq=True)
    
                #Calculate output values: Z, VSWR, Power
                tData.processOutputs(apply_compensation=False)
        
                #Fix a bug from FW assignment of tcpId
                tData.tcpId[tData.tcpId >= 64] -= 64
        
                #Calculate reconstructed scope traces
                tData.process_scope_traces(10)
        
                #Calclate the number of TCPs
                nTcp = np.max(tData.tcpId)+1
        
        
                ##############################
                #Plotting code
                #############################
                if False:
                    import numpy as np
                    import matplotlib.pyplot as plt
                    
                    fig, ax1 = plt.subplots()
                    
                    ax1.plot(tData.scopephi,tData.scopev,color='red',linewidth=2)    
                    #ax1.plot(tData.scopephi,real(mean(tData.Avs[:,0])*exp(1j*tData.scopephi)),color='grey',linewidth=2)    
                    
                    ax1.set_xlabel('Phase (rad)')
                    # Make the y-axis label and tick labels match the line color.
                    ax1.set_ylabel('Voltage (V)', color='red')
                    for tl in ax1.get_yticklabels():
                        tl.set_color('red')
                    
                    
                    ax2 = ax1.twinx()
                    ax2.plot(tData.scopephi,tData.scopei,color='green',linewidth=2)
                    #ax2.plot(tData.scopephi,real(mean(tData.Ais[:,0])*exp(1j*tData.scopephi)),color='grey',linewidth=2)    
                
                    ax2.set_ylabel('Current (A)', color='green')
                    for tl in ax2.get_yticklabels():
                        tl.set_color('green')
                    plt.grid()
                    plt.title("Reconstructed voltage and current waveforms")
                    
                    plt.show()
                    
                    figure()
                    title("First voltage measurement across record number")
                    plot(tData.vmeas[:,0],marker='+')
                    xlabel("Record number")
                    ylabel("Voltage (V)")
                    grid('on')
                
                    figure()
                    title("First current measurement across record number")
                    plot(tData.imeas[:,0],marker='+')
                    xlabel("Record number")
                    ylabel("Current (A)")
                    grid('on')
                    
                    #Plot outputs by TCP
                    newfile=True
                    for i in use_tcps:
                        
                        figure()
                        subplot(211)
                        title("TCP {0}".format(i))
                        plot(np.real(tData.Z[tData.tcpId==i]))
                        ylabel("Zr")
                        grid()
                        subplot(212)
                        plot(np.imag(tData.Z[tData.tcpId==i]))
                        grid()
                        ylabel("Zi")
                        xlabel("Record Number")
                        
                        figure()
                        subplot(211)
                        title("TCP {0}".format(i))
                        plot(abs(tData.V1[tData.tcpId==i]))
                        ylabel("Peak Voltage (V)")
                        grid()
                        subplot(212)
                        plot(abs(tData.I1[tData.tcpId==i]))
                        ylabel("Peak Current (A)")
                        xlabel("Record Number")
                        grid()
                        
                        figure()
                        subplot(311)
                        title("TCP {0}".format(i))
                        plot(tData.PF[tData.tcpId==i])
                        ylabel("Power Forward (W)")
                        grid()
                        
                        subplot(312)        
                        title("TCP {0}".format(i))
                        plot(tData.PR[tData.tcpId==i])
                        ylabel("Power Reflected (W)")
                        xlabel("Record number")
                        grid()
                        show()
                        
                        subplot(313)        
                        title("TCP {0}".format(i))
                        plot(tData.Pload[tData.tcpId==i])
                        ylabel("Load power (W)")
                        xlabel("Record number")
                        grid()
                        show()
