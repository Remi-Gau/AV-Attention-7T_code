#!/bin/bash

# freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0

#export FREESURFER_HOME=/home/SHARED/Program/freesurfer
#export SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data
#source $FREESURFER_HOME/SetUpFreeSurfer.sh


# ?h.G_temp_sup-G_T_transv.label : HG

# ?h.S_temporal_sup.label: STS

# ?h.G_temp_sup-Plan_tempo.label : PT

clear

SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data/

Subjects_List='01 02 03 06 08 09 10 11 12'
Subjects_List='09 11'

HS='lh rh'

for SubjectInd in $Subjects_List
do
	
	echo "\n"
	echo "###################################"
	echo "#	PROCESSING SUBJECT $SubjectInd	#"
	echo "###################################"
	echo "\n"
	
	SourceImage=`ls $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_400MNI.nii`
	
	# coregister .4 mm isometric qT1 image to anatomical surface
	tkregister2 --mov $SourceImage --s $SubjectInd/Structural/FSL_LowDef --regheader --noedit \
	--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat	

	echo "\n\n"
	

	UpSampImage=`ls $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/DownSamp_1mm_UNI_700MNI.nii`
	
	mkdir $SUBJECTS_DIR$SubjectInd/Structural/400MNI/RegData
	
	tkregister2 --mov $UpSampImage --s $SubjectInd/Structural/FSL_LowDef --regheader --noedit \
	--reg $SUBJECTS_DIR$SubjectInd/Structural/400MNI/RegData/register.dat
	
	echo "\n\n"
	
	
	for hs in $HS
	do
	
		# Projects T1 map values at mid cortical depth on surface
		mri_vol2surf --mov $SourceImage --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register.dat --surf white \
		--hemi $hs --projfrac 0.5 --o $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_$hs.mgh
	
		echo "\n\n"
	
	
		# Smooths surface
		mri_surf2surf --s $SubjectInd/Structural/FSL_LowDef --hemi $hs \
		--srcsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_$hs.mgh \
		--trgsurfval $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_${hs}_smooth.mgh \
		--fwhm-trg 4
	
		echo "\n\n"
	
		cp $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/${hs}.G_temp_sup-G_T_transv.label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/HG_${hs}_BU.label
		cp $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/${hs}.G_temp_sup-Plan_tempo.label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/PT_${hs}_BU.label
	
	  	# Visualize to draw ROI # 
	  	tksurfer $SubjectInd/Structural/FSL_LowDef/ $hs inflated -curv \
	  	-patch $hs.temp.patch \
	  	-overlay $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/T1_${hs}_smooth.mgh \
	  	-fminmax 450 900 -foffset 2500 -truncphaseflag 1 -invphaseflag 1
	
		echo "\n\n"

	
		# LABEL 2 VOL
	
		mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/PT_${hs}.label \
		--temp $UpSampImage \
		--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
		--proj frac -.5 1.1 .1 \
		--subject $SubjectInd/Structural/FSL_LowDef \
		--hemi ${hs} \
		--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/ROI/PT_${hs}.nii
		
		echo "\n\n"

		mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/HG_${hs}.label \
		--temp $UpSampImage \
		--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
		--proj frac -.5 1.1 .1 \
		--subject $SubjectInd/Structural/FSL_LowDef \
		--hemi ${hs} \
		--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/ROI/HG_${hs}.nii
	
	done
done
