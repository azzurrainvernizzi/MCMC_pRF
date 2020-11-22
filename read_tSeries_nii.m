function ts = read_tSeries_nii(inputfile)

% where - inputfile is the name of a 3D stimulus in the NifTI format
% returns a file img.mat in the same directory as <inputfile> ready for the
% MCMC pRF computing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if nifti support exists
% otherwise add fastECM-supplied version
if (exist('load_untouch_nii')~=2)
    
    fprintf('Software for reading NifTI images not found.        \n');
    fprintf('Using Tools for Nifti/Analyze for NifTI file I / O. \n');
    fprintf('  www.mathworks.com/matlabcentral/fileexchange/8797 \n');
    usenifti=1;                    
    npath=[fileparts(which('run_MCMC_pRF.m')) filesep 'tools4nifti'];
    addpath(npath);                   
    
else    

  usenifti=0;                      

end 

% load NifTI input file
fprintf('reading %s ...\n',inputfile);
if (inputfile(1)~=filesep)              % if path does not start with a separator
    
    if (    (isunix) | ...              % and -in windows- not with a drive letter or '\\'
            ( (~isunix) & (inputfile(2)~='\') & (inputfile(2)~=':') ) ...
            )                           % make it an absolute path

      inputfile=[pwd filesep inputfile];
    
    end % if isunix
    
end % if inputfile

% get the directory, base name and extension(s)
[fd,fn,fx2]=fileparts(inputfile);       % dir and base file name + extension (may be combined)

% get the data and clear data from file record
M=load_untouch_nii(inputfile);          % read the information from the NifTI file

ts=M.img;

matfile=[fd filesep fn '.mat'];    
save (matfile,'ts');
    
end