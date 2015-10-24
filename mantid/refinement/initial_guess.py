###########################################################################
## Load the simulated S(Q,E) and convolve with experimental resolution  ###
###########################################################################
rootd='/SNSlocal/scratch/jbq/me8t8/data'
LoadNexus(Filename='{0}/elastic.nxs'.format(rootd), OutputWorkspace='elastic')
LoadNexus(Filename='{0}/exp200K.nxs'.format(rootd), OutputWorkspace='exp200K')
parametervalues='0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14'
for K in parametervalues.split():
	LoadSassena(Filename='{0}/K{1}/T_200/fqt_inc_run1_rms2first.h5'.format(rootd,K), TimeUnit=1.0, OutputWorkspace='incSMK{0}'.format(K))
	if K=='0.11':
		Transpose(InputWorkspace='incSMK0.11_fqt.Re', OutputWorkspace='incSMK0.11_fqt.Re')   #from I(Q,t) to I(t,Q)
		Rebin(InputWorkspace='incSMK0.11_fqt.Re', Params=[0.2,0.2,2.0], OutputWorkspace='incSMK0.11_fqt.Re') #Rebin in Q to (0.3, 0.5,..,1.9)
		Transpose(InputWorkspace='incSMK0.11_fqt.Re', OutputWorkspace='incSMK0.11_fqt.Re')  # from I(t,Q) back to I(Q,t)
		Scale(InputWorkspace='incSMK0.11_fqt.Re', factor=0.5, Operation='Multiply', OutputWorkspace='incSMK0.11_fqt.Re')
	SassenaFFT(InputWorkspace='incSMK{0}'.format(K), FFTonlyRealpart=1, DetailedBalance=1, Temp=200)
	Rebin(InputWorkspace='incSMK{0}_sqw'.format(K), Params=[-0.2,0.0004,0.2], OutputWorkspace='incSMK{0}_sqw'.format(K))
	ws=mtd['incSMK{0}_sqw'.format(K)]  # simulated S(Q,E)
	for iQ in range(9):
		Q=0.03+0.02*iQ  #not neccessary but just to remind the Q-value
		SQE=ws.dataY(iQ)
		SQE[499]=SQE[498]   
		SQE[500]=SQE[501]
for K in parametervalues.split():
	ConvolveWorkspaces(Workspace1='elastic', Workspace2='incSMK{0}_sqw'.format(K), OutputWorkspace='simSMK{0}'.format(K))
for K in parametervalues.split():
	Scale(InputWorkspace='simSMK{0}'.format(K), Factor=1.0e-05,Operation='Multiply', OutputWorkspace='simSMK{0}'.format(K))

###########################################################################
#  Use algorith DSFinterp to produce 20 interpolated structure factors ###
###########################################################################
inputworkspaces=['simSMK0.02','simSMK0.03','simSMK0.04','simSMK0.05','simSMK0.06','simSMK0.07',
                 'simSMK0.08','simSMK0.09','simSMK0.10','simSMK0.11','simSMK0.12','simSMK0.13','simSMK0.14',]
parametervalues=[0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14]
targetvalues=[0.03+0.005*i for i in range(20)] #target K varies from 0.03 to 0.13 every 0.005Kcal/mol
outputworkspaces=['interp{0}'.format(i) for i in range(20)] # interpolated structure factors
DSFinterp(ParameterValues=[float(K) for K in parametervalues], Workspaces=inputworkspaces, LoadErrors=0,
          TargetParameters=targetvalues, OutputWorkspaces=outputworkspaces,
          LocalRegression=1, RegressionWindow=4, RegressionType='quadratic')

###########################################################################
# Fit spectrum with index 8 of the simulated structure factors against experimental spectrum 8
###########################################################################
for i in range(20):
        fit_string  ='name=TabulatedFunction,Workspace=elastic,WorkspaceIndex=8,Scaling=1.0,Shift=0.0;'
        fit_string +='name=TabulatedFunction,Workspace=interp{0},WorkspaceIndex=8,Scaling=1.0,Shift=0.0;'.format(i)
        fit_string +='name=LinearBackground,A0=0.0,A1=0.0'
        Fit(fit_string, InputWorkspace='exp200K', WorkspaceIndex=8, StartX=-0.13, EndX=0.10, CreateOutput = 1, Output= 'fit{0}'.format(i))
        ws=mtd['fit{0}_Parameters'.format(i)]
        chi2=ws.row(6)['Value']
        print targetvalues[i], chi2
'name=TabulatedFunction,WorkspaceIndex=0,Scaling=1,Shift=0;name=TabulatedFunction,WorkspaceIndex=0,Scaling=1,Shift=0;name=LinearBackground,A0=0,A1=0'
