#!/bin/bash

# freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0

#export FREESURFER_HOME=/home/SHARED/Program/freesurfer
#export SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

clear

SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data/

Subjects_List='01 02 03 06 08 09 10 11 12'
Subjects_List='03 06 08 09 10 11 12'

for SubjectInd in $Subjects_List
do

	echo "PROCESSING SUBJECT $SubjectInd\n\n"
	
	mkdir $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1
	mkdir $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData
	
	cp $SUBJECTS_DIR$SubjectInd/Structural/400MNI/T1_400MNI.nii $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_400MNI.nii
	cp $SUBJECTS_DIR$SubjectInd/Structural/400MNI/DownSamp_1mm_UNI_700MNI.nii $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/DownSamp_1mm_UNI_700MNI.nii
	#cp $SUBJECTS_DIR$SubjectInd/Structural/S6_MP2RAGE_5_3_TR5000_iPAT=2_T1_Images.nii $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_400MNI_ori.nii

	SourceImage=`ls $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_400MNI.nii`
	UpSampImage=`ls $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/DownSamp_1mm_UNI_700MNI.nii`	
	
	
	echo "\n\n"
	echo '########################################################################'
	echo '#                    REGISTRATION & PROJECTION                         #'
	echo '########################################################################'
	
	# coregister .4 mm isometric qT1 image to anatomical surface
	tkregister2 --mov $SourceImage --s $SubjectInd/Structural/FSL_LowDef --regheader --noedit --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat	
	
	#for i in $(seq 0 .1 1)
  	#do	
		# Projects qT1 values of the middle of the cortex onto surface for lh --fwhm 0.2
		#mri_vol2surf --mov $SourceImage --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat --surf white --hemi lh --projfrac $i --o $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh_$i.mgh
	
		# Smooth surface
		#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi lh \
		#--srcsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh_$i.mgh \
		#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh_${i}_smooth.mgh\
		#--fwhm-trg 2
		
	  	# Projects qT1 values of the middle of the cortex onto surface for rh
	 	#mri_vol2surf --mov $SourceImage --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat --surf white --hemi rh --projfrac $i --o $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh_$i.mgh
	 	
	 	# Smooth surface
		#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi rh \
		#--srcsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh_$i.mgh \
		#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh_${i}_smooth.mgh\
		#--fwhm-trg 2
  	#done
  	
	#mri_vol2surf --mov $SourceImage --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat --surf white --hemi lh --projfrac 0.5 --o $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh.
	
	#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi lh \
	#--srcsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh.mgh \
	#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh_smooth.mgh \
	#--fwhm-trg 2
	
	#mri_vol2surf --mov $SourceImage --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat --surf white --hemi rh --projfrac 0.5 --o $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh.mgh
	
	#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi rh \
	#--srcsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh.mgh \
	#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh_smooth.mgh \
	#--fwhm-trg 2
  	
  	
  	# Smooth curvature and thichness surfaces
  	#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi rh \
	#--srcsurfval thickness --src_type curv \
	#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/surf/thickness_rh_smooth.mgh \
	#--fwhm-trg 2
	
  	#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi lh \
	#--srcsurfval thickness --src_type curv \
	#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/surf/thickness_lh_smooth.mgh \
	#--fwhm-trg 2
	
  	#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi rh \
	#--srcsurfval curv --src_type curv \
	#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/surf/curv_rh_smooth.mgh \
	#--fwhm-trg 2
	
  	#mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi lh \
	#--srcsurfval curv --src_type curv \
	#--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/surf/curv_lh_smooth.mgh \
	#--fwhm-trg 2
	
  	
  	echo "\n\n"
	echo '########################################################################'
	echo '#                       Visualize to draw ROI                          #'
	echo '########################################################################'
  	
  	
  	# Left hemisphere
  	tksurfer $SubjectInd/Structural/FSL_LowDef lh inflated -curv \
  	-overlay $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_lh_0.5.mgh\
  	-truncphaseflag 1 -invphaseflag 1 -foffset 2500 \
  	-fminmax 450 900 -foffset 2500 -truncphaseflag 1 -invphaseflag 1 \
  	
  	
  	#sclv_set_current_threshold_from_percentile min mid max 	
  	#sclv_load_label_value_file fileName field
  	#rotate_brain_x degrees
  	#translate_brain_x distance
  	#labl_load fileName
  	
  	
  	# Right hemisphere  	
  	#tksurfer $SubjectInd/Structural/FSL_LowDef rh inflated \
  	#-overlay $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_rh_0.5.mgh \
  	#-fminmax 450 900 -foffset 2500 -truncphaseflag 1 -invphaseflag 1
  	  	
  	# lh.G_temp_sup-G_T_transv.label : HG
  	# lh.transversetemporal.label: HG and a bit of S
  	
  	# lh.S_temporal_sup.label: STS
  	
  	# lh.G_temp_sup-Plan_tempo.label : PT
  	
  	# lh.S_temporal_transverse.label : S between HG and PT
  	# lh.G_temp_sup-Lateral.label: STG
  	# lh.superiortemporal.label: STG, STS, PT
  	
	
	echo "\n\n"
	echo '########################################################################'
	echo '#                            Create masks                              #'
	echo '########################################################################'
	
	## transform labels into ROIs at 1 mm isometric
	# coregister 1 mm isometric image 
	#tkregister2 --mov $UpSampImage --s $SubjectInd/Structural/FSL_LowDef --regheader --noedit --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat
	
	# transform label into mask
	#mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/qT1/A1c_T1_lh.label \
	#--temp $UpSampImage \
	#--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
	#--fillthresh .5 \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi lh \
	#--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/qT1/A1c_T1_FS_lh.nii
	
	#mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/qT1/A1c_T1_rh.label \
	#--temp $UpSampImage \
	#--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
	#--fillthresh .5 \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi rh \
	#--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/qT1/A1c_T1_FS_rh.nii
	
done
