# !/bin/bash

# freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0

# http://surfer.nmr.mgh.harvard.edu/fswiki/HiResRecon
# for resolution of 0.4 mm, it needs a  version of mri_tessellate that allows more than 10^6 vertices otherwise it will crash

# ssh -X  les-linux-fs3.bham.ac.uk

# export FREESURFER_HOME=~/Programs/freesurfer
# source $FREESURFER_HOME/SetUpFreeSurfer.sh

# TO DO
# 	- FLAGS FOR GPU -use-gpu
#		- http://www.mail-archive.com/freesurfer%40nmr.mgh.harvard.edu/msg22497.html
#	- COMMENTS ON SANITY CHECKS
# 		- http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/OutputData_freeview
#		- http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/TroubleshootingData


clear

echo '\n'
echo '########################################################################'
echo '#                                                                      #'
echo '#                     RUNNING LOW-DEF RECON ALL                        #'
echo '#                                                                      #'
echo '########################################################################'
echo '\n'

#'02 03 04 06 07 08 09 11 12 13 14 15 16'

# Run 03

Subjects_List='03'

Image=UNI.nii

for SubjectInd in $Subjects_List
do
	
	StartTime=`date +%s`

	SUBJECT=FS
	
	SUBJECTS_DIR=/data/AV_Integration_2/Subjects_Data/Subject_$SubjectInd/Structural/

	SourceFolder=/data/AV_Integration_2/Subjects_Data/Subject_$SubjectInd/Structural/

	## First: downsample high resolution data to 1 mm and process completely

	echo '\n'
	echo '#################################################################'
	echo '#                                                               #'
	echo "#       Subject $SubjectInd : RUNNING LOW-DEF RECON ALL         #"
	echo '#                                                               #'
	echo '#################################################################'
	echo '\n'

	echo "\n\nProcessing subject $SubjectInd \n\n"
		 
	recon-all -motioncor -talairach -tal-check -i $SourceFolder$Image -s $SUBJECTS_DIR$SUBJECT
	# If the talairach registration fails an automated correction is being performed but if for whatever reason the result is not being used, try:
        #cp $SourceFolder/FS/mri/transforms/talairach.auto.xfm \
        # $SourceFolder/FS/mri/transforms/talairach.xfm
	# recon-all -tal-check -s $SUBJECT 
        # If it still fails try looking at the manual talairach registration in the wiki of Freesurfer.

	mri_nu_correct.mni \
	--i $SUBJECTS_DIR$SUBJECT/mri/orig.mgz \
	--o $SUBJECTS_DIR$SUBJECT/mri/nu.mgz \
	--proto-iters 1000 --distance 15 --fwhm 0.15 --n 1 --uchar \
	$SUBJECTS_DIR$SUBJECT/mri/transforms/talairach.xfm

	recon-all -mprage -normalization -skullstrip -s $SUBJECTS_DIR/$SUBJECT

	########################################################################
	#                                                                      #
	#            Visualize brainmask.mgz before continuing	               #
	#                                                                      #
	########################################################################
 

	# freeview -v $SourceFolder/FS/mri/orig/001.mgz $SourceFolder/FS/mri/brainmask.mgz
	
	# If skullstripping does not look satisfactory, try:
	# mri_watershed -T1 -atlas -h 35 -brain_atlas $FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca $SourceFolder/FS/mri/transforms/talairach_with_skull.lta SourceFolder/FS/mri/T1.mgz SourceFolder/FS/mri/brainmask.auto.mgz 
	
	# If small structures are still attached to the brain, try:
	# mri_gcut -110 -mult $SourceFolder/$SourceFolder/FS/mri/brainmask.auto.mgz $SourceFolder/$SourceFolder/FS/mri/T1.mgz $SourceFolder/FS/mri/brainmask.auto.mgz 

	# More information: http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/SkullStripFix_freeview


	# If skullstriping was satisfactory: 
	cp $SUBJECTS_DIR$SUBJECT/mri/brainmask.auto.mgz $SUBJECTS_DIR$SUBJECT/mri/brainmask.mgz
	recon-all -autorecon2 -autorecon3 -mprage -s $SUBJECTS_DIR$SUBJECT


	# Mail results
	tail -4 $SUBJECTS_DIR$SUBJECT/scripts/recon-all.log | mailx -s "Analysis LD of subject $SubjectInd finished" -c remi.gau@gmail.com remi_gau@hotmail.com	


	EndTime=`date +%s`
	echo execution time was `expr $EndTime - $StartTime` s.

	
done	
