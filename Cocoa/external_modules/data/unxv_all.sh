#!/bin/bash

rm -rf ./sn_data/
tar xf sn_data.xz

rm -rf ./des_data/
tar xf des_data.xz

rm -rf ./bicep_keck_2015/
tar xf bicep_keck_2015.xz

rm -rf ./bao_data/
tar xf bao_data.xz

rm -rf ./simons_observatory/
tar xf simons_observatory.xz

rm -rf ./spt_hiell_2020/
tar xf spt_hiell_2020.xz

rm -rf ./h0licow_distance_chains/
tar xf h0licow_distance_chains.xz

if [ -z "${THREAD_UNXZ}" ]; then
	# ---------------------------------------------
	cd ./planck
	# ---------------------------------------------
	rm -rf ./spt3g_Y1_EETE.clik/
	tar xf spt3g_Y1_EETE.clik.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING CAMSPEC (2018)'
	rm -rf ./CamSpec2018/
	tar xf CamSpec2018.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES'
	rm -rf ./planck_supp_data_and_covmats/
	tar xf planck_supp_data_and_covmats.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0)'
	rm -rf ./plc_2.0/
	tar xf plc_20.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0)'
	cd ./plc_3.0
	rm -rf ./lensing/
	tar xf lensing.xz
	rm -rf ./low_l/
	tar xf low_l.xz
	cd ./hi_l
	rm -rf ./plik/
	tar xf plik.xz
	rm -rf ./plik_lite/
	tar xf plik_lite.xz
	# ---------------------------------------------
	cd ./camspec
	rm -rf ./camspec_10.7HM_1400_TT_small.clik/
	tar xf camspec_107HM_1400_TT_smallclik.xz
	rm -rf ./camspec_10.7HM_1400_TTTEEE.clik/
	tar xf camspec_10.7HM_1400_TTTEEEclik.xz
	# ---------------------------------------------
else
	# ---------------------------------------------
	cd ./planck
	# ---------------------------------------------
	rm -rf ./spt3g_Y1_EETE.clik/
	tar xf spt3g_Y1_EETE.clik.xz &
	proc10=$!
	# ---------------------------------------------
	rm -rf ./CamSpec2018/
	tar xf CamSpec2018.xz &
	# ---------------------------------------------
	proc1=$!
	rm -rf ./planck_supp_data_and_covmats/
	tar xf planck_supp_data_and_covmats.xz &
	proc2=$!
	# ---------------------------------------------
	rm -rf ./plc_2.0/
	tar xf plc_20.xz &
	proc3=$!
	# ---------------------------------------------
	cd ./plc_3.0
	rm -rf ./lensing/
	tar xf lensing.xz &
	proc4=$!
	rm -rf ./low_l/
	tar xf low_l.xz &
	proc5=$!
	cd ./hi_l
	rm -rf ./plik/
	tar xf plik.xz &
	proc6=$!
	rm -rf ./plik_lite/
	tar xf plik_lite.xz	&
	proc7=$!
	# ---------------------------------------------
	cd ./camspec
	rm -rf ./camspec_10.7HM_1400_TT_small.clik/
	tar xf camspec_107HM_1400_TT_smallclik.xz &
	proc8=$!
	rm -rf ./camspec_10.7HM_1400_TTTEEE.clik/
	tar xf camspec_10.7HM_1400_TTTEEEclik.xz &
	proc9=$!
	# ---------------------------------------------
	# ---------------------------------------------
	# ---------------------------------------------
	echo 'DECOMPRESSING CAMSPEC (2018)'
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES'
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0)'
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0)'
	echo 'DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH'
	wait "$proc1" "$proc2" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9" "$proc10"
fi
