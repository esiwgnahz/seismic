include ${SEPINC}/SEP.top

include EXEMakefile
include VOLVEMakefile
wavelet.H:
	Wavelet wavelet=ricker2 fund=25 tdelay=0.13 n1=1000 d1=0.001 > $@

wavelet3D.H:
	Wavelet wavelet=ricker2 fund=25 tdelay=0.13 n1=1000 d1=0.0009 > $@


3Dvel.H:
	Spike n1=50 n2=201 n3=201 | Add scale=2800 > 1.H
	Spike n1=50 n2=201 n3=201 | Add scale=3200 > 2.H
	Spike n1=50 n2=201 n3=201 | Add scale=3600 > 3.H
	Cat axis=1 1.H 2.H 3.H > $@
	echo d1=10 d2=10 d3=10 >> $@
	Window3d n3=1 < $@ > 2Dvel.H

3Dref.H: 3Dvel.H
	Window3d f1=1 < 3Dvel.H | Pad end1=1 extend=1 > junk.H
	Math file1=junk.H file2=3Dvel.H exp="(file1-file2)/800" > $@
	echo o1=0. >> $@

3Drho.H:
	Spike n1=50 n2=201 n3=201 | Add scale=1 > 1.H
	Spike n1=50 n2=201 n3=201 | Add scale=1.2 > 2.H
	Spike n1=50 n2=201 n3=201 | Add scale=1.5 > 3.H
	Cat axis=1 1.H 2.H 3.H > $@
	echo d1=10 d2=10 d3=10 >> $@
	Window3d n3=1 < $@ > 2Drho.H

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

rhodata1%.H: wavelet.H 3Dvel.H 3Drho.H 
	${MOD} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=1 withRho=1 rho=3Drho.H ${param$*} fmax=60 data=$@ 

3Drhodata1%.H: wavelet.H 3Dvel.H 3Drho.H 
	${MOD} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=0 withRho=1 rho=3Drho.H ${param$*} fmax=60 data=$@ 

3Ddata1%.H: wavelet3D.H 3Dvel.H 3Drho.H 
	${MOD} stdin=wavelet3D.H vel=3Dvel.H wave=wave1.H twoD=0 withRho=0 rho=3Drho.H ${param$*} fmax=60 data=$@ 

3DBorndata1%.H: wavelet3D.H 3Dvel.H 3Drho.H 3Dref.H
	${BMOD} stdin=wavelet3D.H vel=3Dvel.H wave=wave1.H twoD=0 ${param$*} fmax=60 data=$@ ref=3Dref.H

3Drtm1%.H: wavelet3D.H 3Dvel.H 3Drho.H 3Dref.H
	${RTM} stdin=wavelet3D.H vel=3Dvel.H wave=wave1.H twoD=0 ${param$*} fmax=60 data=$@ ref=3Dref.H

Borndata1%.H: wavelet.H 3Dvel.H 3Drho.H 3Dref.H
	${BMOD} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=1 ${param$*} fmax=60 data=$@ ref=3Dref.H

rhodata2%.H: traces.H 3Drho.H
	${REF} \
	land=1 \
	vel=2Dvel.H \
	dens=2Drho.H \
	vel_dens_mask=2Dvel.H \
	traces=traces.H \
	density=1 \
	fmax=60  \
	task=1 \
	${param$*} \
	stdin=wavelet.H > $@

BORNdata2%.H: traces.H 3Drho.H
	Window3d n3=1 < 3Dref.H > 2Dref.H
	echo d1=10 d2=10 >> 2Dref.H
	${BORNREF} \
	land=1 \
	Reflectivity=2Dref.H \
	vel=2Dvel.H \
	vel_dens_mask=2Dvel.H \
	traces=traces.H \
	fmax=60  \
	snap_coef=0.5 \
	task=1 \
	${param$*} tmpwave1=tmpwave1.H tmpwave2=tmpwave2.H \
	stdin=wavelet.H Modeling=$@ 

rhodata3%.H:  3Drho.H
	${BRU} \
	stdin=wavelet.H \
	vel=2Dvel.H \
	rho=2Drho.H \
	fmax=60 \
	max_file_size=10 \
	shot_z=0 \
	ishot_x=101 \
	ishot_y=100 \
	${param$*} \
	Seismogram=junk1.H
	Transp plane=12 < junk1.H > $@

include ${SEPINC}/SEP.bottom
