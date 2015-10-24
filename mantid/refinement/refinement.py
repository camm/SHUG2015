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
guess={'Scaling':1.0, 'Intensity':1.0, 'K':0.11}
inputworkspaces=' '.join( ['simSMK{0}'.format(K) for K in parametervalues.split()])
for iw in range(8,-1,-1): 
        Q=0.3+iw*0.2
	fit_string = 'name=TabulatedFunction,Workspace=elastic,WorkspaceIndex={0},Scaling={1},constraints=(0.0001<Scaling);'.format(iw,guess['Scaling'])+\
		'name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(inputworkspaces,parametervalues) +\
		'LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=4,' +\
		'WorkspaceIndex={0},Intensity={1},TargetParameter={2},'.format(iw,guess['Intensity'],guess['K']) +\
		'constraints=(0.0001<Intensity);' +\
		'name=LinearBackground,A0=0.0,A1=0.0'
	Fit(fit_string, InputWorkspace='exp200K', WorkspaceIndex=iw, StartX=-0.13, EndX=0.1, CreateOutput = 1, Output='fit{0}'.format(iw) )
        ws=mtd['fit{0}_Parameters'.format(iw)]
        guess['Scaling']=ws.row(0)['Value']
        guess['Intensity']=ws.row(2)['Value']
        Koptimal=ws.row(3)['Value']
        print 'Q=', Q, 'L={0:5.1f}'.format( 6.2832/Q), 'K={0:6.4f}'.format(Koptimal), 'Chi2={0:3.1f}'.format(ws.row(6)['Value'])
