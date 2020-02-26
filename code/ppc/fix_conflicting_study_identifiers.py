#go thru all dicoms and look for inconsistent studyUIDs, if there are any, change them to match the first one.
#run this on each subject's dicom directory that had an interrupted scanning session
import pydicom
import os
import glob
tmpdir=os.getcwd()
alldcm = glob.glob(tmpdir + '/ses-baselineYear1Arm1/*/*/*.dcm')
for jj in range(0,len(alldcm)):
    ds = pydicom.dcmread(alldcm[jj])
    if jj is 0:
        studyUID = ds.StudyInstanceUID
    ds.StudyInstanceUID= studyUID
    ds.save_as(alldcm[jj])
