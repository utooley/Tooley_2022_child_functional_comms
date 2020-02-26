import sys
import json
import bisect
from glob import glob
from os.path import join, splitext
from bids.grabbids import BIDSLayout
from dateutil.parser import parse

# subj_dir *must* have trailing /
subj_dir = '/data/jux/mackey_group/public_data/ABCD/release2_site20/niftis/'
subj=sys.argv[1]
#subj=
sess = '01'
data_suffix = '.nii.gz'

layout = BIDSLayout(subj_dir)


def files_to_dict(file_list):
    """Convert list of BIDS Files to dictionary where key is
    acquisition time (datetime.datetime object) and value is
    the File object.
    """
    out_dict = {}
    for f in file_list:
        fn = f.filename
        with open(fn, 'r') as fi:
            data = json.load(fi)
        dt = parse(data['AcquisitionTime'])
        out_dict[dt] = f
    return out_dict

# Get json files for field maps
fmap_jsons = layout.get(subject= subj, modality='fmap', extensions='.json')

for dir_ in ['AP', 'PA']:
    # Run field map directions independently
    dir_jsons = [fm for fm in fmap_jsons if '_dir-{0}_'.format(dir_) in fm.filename]
    fmap_dict = files_to_dict(dir_jsons)
    dts = sorted(fmap_dict.keys())

    intendedfor_dict = {fmap.filename: [] for fmap in dir_jsons}

    # Get all scans with associated field maps (bold + dwi)
    func_jsons = layout.get(subject= subj, type='bold', extensions='.json') +\
                 layout.get(subject= subj, type='dwi', extensions='.json')
    func_dict = files_to_dict(func_jsons)
    for func in func_dict.keys():
        fn, _ = splitext(func_dict[func].filename)
        fn += data_suffix
        fn=fn.split(subj_dir)[-1]
        fn = fn.split("/",1)[1]
        # Find most immediate field map before scan
        idx = bisect.bisect_right(dts, func) - 1
        fmap_file = fmap_dict[dts[idx]].filename
        intendedfor_dict[fmap_file].append(fn) #should this be adjusted so that it's immediately before or after?

    for fmap_file in intendedfor_dict.keys():
        with open(fmap_file, 'r') as fi:
            data = json.load(fi)

        # No overwriting, for now
        if 'IntendedFor' not in data.keys():
            data['IntendedFor'] = intendedfor_dict[fmap_file]
            with open(fmap_file, 'w') as fo:
                json.dump(data, fo)
