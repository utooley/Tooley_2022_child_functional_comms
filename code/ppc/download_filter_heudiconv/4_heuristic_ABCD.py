import os

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes

def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where

    allowed template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """

    # paths done in BIDS format
    t1w = create_key('sub-{subject}/anat/sub-{subject}_run-{item:02d}_T1w')
    t2w = create_key('sub-subject}/anat/sub-{subject}_run-{item:02d}_T2w')
    rest = create_key('sub-{subject}/func/sub-{subject}_task-rest_run-{item:02d}_bold')
    #task_nback = create_key('sub-{subject}/func/sub-{subject}_task-nback_run-{item:02d}_bold')
    #task_num = create_key('sub-{subject}/func/sub-{subject}_task-number_run-{item:02d}_bold')
    dwi = create_key('sub-{subject}/dwi/sub-{subject}_run-{item:02d}_dwi')
    fmap = create_key('sub-{subject}/fmap/sub-{subject}_acq-{label}_dir-{dir}_run-{item:02d}_epi')

    info = {t1w: [], t2w: [], rest: [], dwi:[], fmap: []}

    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:

        * total_files_till_now
        * example_dcm_file
        * series_number
        * dcm_dir_name
        * unspecified2
        * unspecified3
        * dim1
        * dim2
        * dim3
        * dim4
        * TR
        * TE
        * protocol_name
        * is_motion_corrected
        * is_derived
        * patient_id
        * study_description
        * referring_physician_name
        * series_description
        * image_type
        """
        if 'MPR' in s.protocol_name and 'NORM' in s.series_description:
            info[t1w].append(s.series_id)
        if 'T2' in s.protocol_name:
            info[t2w].append(s.series_id)
        if 'rest' in s.protocol_name:
            info[rest].append(s.series_id)
        if 'dMRI' in s.protocol_name and (s.dim4==103):
            info[dwi].append(s.series_id)
        if 'fMRI' in s.protocol_name and'DistortionMap' in s.protocol_name and '_PA' in s.protocol_name:
            #info[fmap_dwi].append({'item': s.series_id, 'dir': 'AP', 'acq': acq})
            info[fmap].append({'item': s.series_id, 'dir': 'PA', 'label': 'fMRI'});
        if 'fMRI' in s.protocol_name and 'DistortionMap' in s.protocol_name and '_AP' in s.protocol_name:
            info[fmap].append({'item': s.series_id, 'dir': 'AP', 'label': 'fMRI'});
        if 'dMRI' in s.protocol_name and 'DistortionMap' in s.protocol_name and '_AP' in s.protocol_name:
            info[fmap].append({'item': s.series_id, 'dir': 'AP', 'label': 'dMRI'});
        if 'dMRI' in s.protocol_name and 'DistortionMap' in s.protocol_name and '_PA' in s.protocol_name:
            info[fmap].append({'item': s.series_id, 'dir': 'PA', 'label': 'dMRI'});

    return info
