#!/bin/bash

Subjects_List='02 03 04 06 07 08 09 11 12 13 14 15 16'

Folder_List='Analysis Behavioral FieldMap Nifti Raw ROI_MIPAV ROI_MNI Structural'

echo "Creating Logs"

DATETAG="$(date +%Y_%m_%d_%H_%M)"

BackUpLog="/media/BackUp/AV_Integration_7T_2/"$DATETAG"_BackUpLog.txt"
ErrorLog="/media/BackUp/AV_Integration_7T_2/"$DATETAG"_BackUpErrorLog.txt"

touch $BackUpLog
touch $ErrorLog

sleep 5

xterm -hold -e "tail -f ${ErrorLog}" & xterm -hold -e "tail -f ${BackUpLog}" &


echo "Backing up scripts\n"
rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/*.m /media/BackUp/AV_Integration_7T_2/ >>$BackUpLog 2>>$ErrorLog
rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/*.sh /media/BackUp/AV_Integration_7T_2/ >>$BackUpLog 2>>$ErrorLog
rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/*.LayoutXML /media/BackUp/AV_Integration_7T_2/ >>$BackUpLog 2>>$ErrorLog
rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/*.ods /media/BackUp/AV_Integration_7T_2/ >>$BackUpLog 2>>$ErrorLog
rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/*.mat /media/BackUp/AV_Integration_7T_2/ >>$BackUpLog 2>>$ErrorLog

if [ ! -d "/media/BackUp/AV_Integration_7T_2/SubFun" ]
then
	mkdir /media/BackUp/AV_Integration_7T_2/SubFun
fi
rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/SubFun/ /media/BackUp/AV_Integration_7T_2/SubFun/ >>$BackUpLog 2>>$ErrorLog

if [ ! -d "/media/BackUp/AV_Integration_7T_2/Subjects_Data" ]
then
	mkdir /media/BackUp/AV_Integration_7T_2/Subjects_Data/
fi

for SubjectInd in $Subjects_List
do

	echo "\nBacking up subject $SubjectInd"
	
	if [ ! -d "/media/BackUp/AV_Integration_7T_2/Subjects_Data/Subject_$SubjectInd" ]
	then
		mkdir /media/BackUp/AV_Integration_7T_2/Subjects_Data/Subject_$SubjectInd
	fi
	
		rsync -rtvhsc --progress --delete --delete-excluded --itemize-changes /home/SHARED/Experiment/AV_Integration_7T_2/Subjects_Data/Subject_$SubjectInd/ /media/BackUp/AV_Integration_7T_2/Subjects_Data/Subject_$SubjectInd/ >>$BackUpLog 2>>$ErrorLog


done




