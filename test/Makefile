MOD=../bin/LINUX/Acoustic_modeling.x
REF=../../FWI/Progs/GIS/trunk/WI/build/Release/WIA/WI/WIA.x
BRU=/work/FWI/Progs/GIS/trunk/RTM/Acoustic/bin/Acoustic3d.main.x

wavelet.H:
	Wavelet wavelet=ricker2 fund=25 tdelay=0.13 n1=1000 d1=0.001 > $@


3Dvel.H:
	Spike n1=50 n2=201 n3=201 | Add scale=2800 > 1.H
	Spike n1=50 n2=201 n3=201 | Add scale=3200 > 2.H
	Spike n1=50 n2=201 n3=201 | Add scale=3600 > 3.H
	Cat axis=1 1.H 2.H 3.H > $@
	echo d1=10 d2=10 d3=10 >> $@
	Window3d n3=1 < $@ > 2Dvel.H

traces.H:
	echo n1=1000 o1=0. d1=0.001 >> $@
	echo n2=201 o2=0. d2=10 >> $@
	echo n3=1   o3=1000. d3=10 >> $@

param1= surf_type=1 shot_type=0 rec_type=0
param0= surf_type=0 shot_type=0 rec_type=0
param2= surf_type=0 shot_type=0 rec_type=1
param3= surf_type=0 shot_type=1 rec_type=1
param4= surf_type=0 shot_type=1 rec_type=0

data1%.H: wavelet.H 3Dvel.H
	${MOD} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=1 withRho=0 ${param$*} fmax=60 data=$@ 


data2%.H: traces.H
	${REF} \
	land=1 \
	vel=2Dvel.H \
	traces=traces.H \
	fmax=60  \
	task=1 \
	${param$*} \
	stdin=wavelet.H > $@

data3%.H: 
	${BRU} \
	stdin=wavelet.H \
	vel=2Dvel.H \
	fmax=60 \
	max_file_size=10 \
	shot_z=0 \
	ishot_x=101 \
	ishot_y=100 \
	${param$*} \
	Seismogram=junk1.H
	Transp plane=12 < junk1.H > $@
