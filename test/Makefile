include ${SEPINC}/SEP.top

RESDIR=Fig
include EXEMakefile
include VOLVEMakefile
include MARMMakefile
include BPMakefile
include 1DMakefile
include FIGMakefile
#include MATCHMakefile

OPSTEXFLAGS= color=y fat=1 fatmult=1.5 invras=n force= 

wavelet.H:
	Wavelet wavelet=ricker2 fund=25 tdelay=0.13 n1=1000 d1=0.001 > $@

waveletlong.H:
	Wavelet wavelet=ricker2 fund=25 tdelay=0.13 n1=5000 d1=0.001 > $@

wavelet3D.H:
	Wavelet wavelet=ricker2 fund=25 tdelay=0.13 n1=1000 d1=0.0009 > $@


3Dvel.H:
	Spike n1=50 n2=201 n3=21 | Add scale=2800 > 1.H
	Spike n1=50 n2=201 n3=21 | Add scale=3200 > 2.H
	Spike n1=50 n2=201 n3=21 | Add scale=3600 > 3.H
	Cat axis=1 1.H 2.H 3.H > $@
	echo d1=10 d2=10 d3=10 >> $@
	echo o1=0. o2=0. o3=-100 >> $@
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

1shot3D.H:
	Spike n1=1000 n2=4221 >> $@
	echo d1=0.001 o1=0 >> $@
	echo d2=10 o2=0. >> $@
	echo d3=10 o3=0. >> $@

1shotheader3D.H:
	Spike n1=1 n2=201 n3=21 > 1.H
	Math file1=1.H exp="0." > 0.H
	Math file1=1.H exp="1000." > sx.H
	echo d1=1 d2=10 d3=10 o1=0. o2=0. o3=0. >> 1.H
	Transp plane=12 reshape=1,3 < 1.H | Tpow tpow=1 > 2.H
	Math file1=2.H exp="file1-10" | Transp plane=12 reshape=2,3> rx.H
	Transp plane=13 < 1.H | Tpow tpow=1 > 2a.H
	Math file1=2a.H exp="file1-110" | Transp plane=13 > ry.H
	echo n1=1 n2=4221 n3=1 >> rx.H
	echo n1=1 n2=4221 n3=1 >> ry.H
	echo n1=1 n2=4221 n3=1 >> sx.H
	echo n1=1 n2=4221 n3=1 >> 0.H
	Cat axis=1 sx.H 0.H 0.H rx.H ry.H 0.H > $@

5shots.H:
	Spike n1=5000 n2=2700 n3=1 >> $@
	echo d1=0.001 o1=0 >> $@
	echo d2=10 o2=0. >> $@
	echo d3=10 o3=0. >> $@

5shotsheader.H:
	Spike n1=1 n2=400 > 1.H
	Math file1=1.H exp="0." > 0.H
	echo d1=1 d2=10 o1=0. o2=0. >> 1.H
	Transp plane=12 < 1.H | Tpow tpow=1 > 2.H
	Math file1=2.H exp="file1-10" | Transp plane=12 > rx.H
	Math file1=1.H exp="1000." > sx0.H
	Math file1=1.H exp="1120." > sx1.H
	Math file1=1.H exp="1210." > sx2.H
	Math file1=1.H exp="1350." > sx3.H
	Math file1=1.H exp="1467." > sx4.H
	Cat axis=1 sx0.H 0.H 0.H rx.H 0.H 0.H > junk0.H
	Cat axis=1 sx1.H 0.H 0.H rx.H 0.H 0.H > junk1.H
	Cat axis=1 sx2.H 0.H 0.H rx.H 0.H 0.H > junk2.H
	Cat axis=1 sx3.H 0.H 0.H rx.H 0.H 0.H > junk3.H
	Cat axis=1 sx4.H 0.H 0.H rx.H 0.H 0.H > junk4.H
	Cat axis=2 junk0.H junk1.H junk2.H junk3.H junk4.H > junk5.H
	Spike n1=1 n2=700 > 1.H
	Math file1=1.H exp="0." > 0.H
	echo d1=1 d2=10 o1=0. o2=0. >> 1.H
	Transp plane=12 < 1.H | Tpow tpow=1 > 2.H
	Math file1=2.H exp="file1-10" | Transp plane=12 > rx.H
	Math file1=1.H exp="1800." > sx0.H
	Cat axis=1 sx0.H 0.H 0.H rx.H 0.H 0.H > junk6.H
	Cat axis=2 junk6.H junk5.H > $@

1shot.H:
	Spike n1=1000 n2=201 n3=1 >> $@
	echo d1=0.001 o1=0 >> $@
	echo d2=10 o2=0. >> $@
	echo d3=10 o3=0. >> $@

1shotheader.H:
	Spike n1=1 n2=201 > 1.H
	Math file1=1.H exp="0." > 0.H
	Math file1=1.H exp="1000." > sx.H
	echo d1=1 d2=10 o1=0. o2=0. >> 1.H
	Transp plane=12 < 1.H | Tpow tpow=1 > 2.H
	Math file1=2.H exp="file1-10" | Transp plane=12 > rx.H
	Cat axis=1 sx.H 0.H 0.H rx.H 0.H 0.H > $@

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

3DBorndata.H: wavelet.H 3Dvel.H 3Dref.H 1shotheader3D.H 1shot3D.H
	${BMODSEP} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=0 Born=1 fmax=60 data=$@ ref=3Dref.H traces=1shot3D.H coordfile=1shotheader3D.H fmax=16 aperture_x=2 aperture_y=200 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=6 \
	withRho=0 \
	snapi=1 ntaper=20 num_threads=8 \
	data=$@ maxsize=8

Borndata.H: wavelet.H 3Dvel.H 3Dref.H 1shotheader.H 1shot.H
	Window3d n3=1 < 3Dref.H > 2Dref.H
	Window3d n3=1 < 3Dvel.H > 2Dvel.H
	${BMODSEP} stdin=wavelet.H vel=2Dvel.H wave=wave1.H twoD=1 Born=1 fmax=60 data=$@ ref=2Dref.H traces=1shot.H coordfile=1shotheader.H fmax=16 aperture_x=2 aperture_y=1 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=6 \
	withRho=0 \
	snapi=1 ntaper=20 num_threads=8 \
	data=$@ maxsize=1 

3DMIGBorndata.H: 3DBorndata.H wavelet.H 3Dvel.H 3Dref.H 1shotheader.H
	${LSRTMSEP} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=0 Born=1 fmax=60 data=$@ ref=3Dref.H traces=3DBorndata.H coordfile=1shotheader3D.H fmax=16 aperture_x=2 aperture_y=1 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	snapi=1 ntaper=20 num_threads=8 \
	image=$@ maxsize=1 LSRTM=1 

MIGBorndata.H: Borndata.H wavelet.H 3Dvel.H 3Dref.H 1shotheader.H
	Window3d n3=1 < 3Dref.H > 2Dref.H
	Window3d n3=1 < 3Dvel.H > 2Dvel.H
	${LSRTMSEP} stdin=wavelet.H vel=2Dvel.H wave=wave1.H twoD=1 Born=1 fmax=60 data=$@ ref=2Dref.H traces=Borndata.H coordfile=1shotheader.H fmax=16 aperture_x=2 aperture_y=1 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	snapi=1 ntaper=20 num_threads=8 \
	image=$@ maxsize=1 LSRTM=1 

3DBornrandom.H: wavelet.H 3Dvel.H 3Dref.H 1shotheader.H 1shot.H
	< 3Dref.H  sfnoise rep=1 > 3Drefnoise.H
	${BMODSEP} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=0 Born=1 fmax=60 data=$@ ref=3Drefnoise.H traces=1shot.H coordfile=1shotheader.H fmax=16 aperture_x=2 aperture_y=200 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	snapi=4 ntaper=20 num_threads=8 \
	data=$@ maxsize=10 

3DMIGBornrandom.H: 3DBornrandom.H wavelet.H 3Dvel.H 3Dref.H 1shotheader.H
	sfnoise rep=1 < 3DBornrandom.H > noisedat.H
	${LSRTMSEP} stdin=wavelet.H vel=3Dvel.H wave=wave1.H twoD=0 Born=1 fmax=60 data=$@ ref=3Dref.H traces=noisedat.H coordfile=1shotheader.H fmax=16 aperture_x=2 aperture_y=200 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	snapi=4 ntaper=20 num_threads=8 \
	image=img.H maxsize=10 LSRTM=1 
	Window3d n4=1 < img.H > $@

Bornrandom.H: wavelet.H 3Dvel.H 3Dref.H 1shotheader.H 1shot.H
	Window3d n3=1 < 3Dref.H | Noise rep=1 > 2Drefnoise.H
	Window3d n3=1 < 3Dvel.H > 2Dvel.H
	${BMODSEP} stdin=wavelet.H vel=2Dvel.H wave=wave1.H twoD=1 Born=1 fmax=60 data=$@ ref=2Drefnoise.H traces=1shot.H coordfile=1shotheader.H fmax=16 aperture_x=2 aperture_y=1 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	snapi=40 ntaper=20 num_threads=8 \
	data=$@ maxsize=1 

MIGBornrandom.H: Borndata.H wavelet.H 3Dvel.H 3Dref.H 1shotheader.H
	Window3d n3=1 < 3Dref.H > 2Dref.H
	Window3d n3=1 < 3Dvel.H > 2Dvel.H
	Noise rep=1 < Borndata.H > noisedat.H
	${LSRTMSEP} stdin=wavelet.H vel=2Dvel.H wave=wave1.H twoD=1 Born=1 fmax=60 data=$@ ref=2Dref.H traces=noisedat.H coordfile=1shotheader.H fmax=16 aperture_x=2 aperture_y=1 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	snapi=40 ntaper=20 num_threads=8 \
	image=img.H maxsize=1 LSRTM=1 
	Window3d n3=1 < img.H > $@

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
	${BORNREF} \
	land=1 \
	Reflectivity=2Dref.H \
	vel=2Dvel.H \
	vel_dens_mask=2Dvel.H \
	traces=traces.H \
	fmax=16  \
	snap_coef=0. \
	task=1 \
	${param$*} tmpwave1=tmpwave1.H tmpwave2=tmpwave2.H \
	stdin=wavelet.H Modeling=$@ 

MIGBORNdata2%.H: BORNdata2%.H traces.H 3Drho.H
	${BORNREF} \
	land=1 \
	Reflectivity=2Dref.H \
	vel=2Dvel.H \
	vel_dens_mask=2Dvel.H \
	traces=BORNdata2$*.H \
	fmax=16  \
	snap_coef=2.5 dottest=1 \
	task=1 \
	${param$*} tmpwave1=tmpwave1.H tmpwave2=tmpwave2.H \
	stdin=wavelet.H RTM=$@ 

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

testgeo%.H: 5shotsheader.H 5shots.H waveletlong.H 2Dvel.H
	time ${AFWI} stdin=waveletlong.H vel=2Dvel.H twoD=1 \
	fmax=60 \
	traces=5shots.H \
	coordfile=5shotsheader.H \
	fmax=16 aperture_x=2 aperture_y=1 \
	keygx=4 \
	keygy=5 \
	keysx=1 \
	keysy=2 \
	keyselev=3 \
	keygelev=3 \
	withRho=0 \
	grad=grad$*.H \
	illu=illu$*.H \
	snapi=$* ntaper=20 num_threads=6 \
	data=$@ maxsize=1 > /dev/null

include ${SEPINC}/SEP.bottom
