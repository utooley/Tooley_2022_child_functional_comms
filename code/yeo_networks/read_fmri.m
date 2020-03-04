function vol = read_fmri(fmri_name)

% [fmri, vol] = read_fmri(fmri_name)
% Given the name of averaged correlation profile file (fmri_name), this
% function read in the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%       It can be .nii.gz file for correlation profile in fsaverage*
%       spaces; or .mat file for correlation profile in fs_LR_32k space
%       (there must be a variable called 'profile_mat' in the .mat file).
%
% Output:
%     - vol:
%       A num_voxels x num_timepoints matrix which is the reshaped 'vol'
%       structure for .nii.gz file or the variable 'profile_mat' for .mat file.
%

if (~isempty(strfind(fmri_name, '.nii.gz')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = fmri.vol;
    vol_size = size(vol);
    if(length(vol_size) < 4)
        vol = reshape(vol, prod(vol_size(1:3)), 1);
    else
        vol = reshape(vol, prod(vol_size(1:3)), vol_size(4));
    end
    fmri.vol = [];
else
    load(fmri_name);
    vol = profile_mat;
end
