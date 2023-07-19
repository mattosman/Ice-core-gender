# Ice core authorship gender analysis (ca. 1969-2021 CE)

This repository contains the MATLAB scripts and underlying source data needed to reproduce the ice core science authorship gender analyses described in Koffman *et al.* (*in revision*). <br>

To generate these analyses, users must first run “saveAbstracts.m” and  “savePersons.m” which download, process, organize, then save five decades worth of ice core-related abstract meta-data into a more convenient MATLAB-based data format.  The generated data files, “abstracts.mat” and “persons.mat”, are each approximately 75 mB.  Both .mat files are needed to run the various plotting scripts located in the subdirectory “Figures”. <br>

To reproduce Figures 1-3 and Extended Data Figures (EDF) 1-2 of Koffman et al. (in revision):
```matlab
% make sure “Ice-core-gender” is your home directory …
saveAbstracts % creates abstracts.mat
savePersons % creates persons.mat
cd Figures
	Fig1 % reproduces Koffman et al. (in revision) Fig. 1
	Fig2 % reproduces Koffman et al. (in revision) Fig. 2
	Fig3 % reproduces Koffman et al. (in revision) Fig. 3
	EDF1 % reproduces Koffman et al. (in revision) EDF 1
	EDF2 % reproduces Koffman et al. (in revision) EDF 2
cd ../
```

Please note this repository includes the subdirectories “cbrewer” forked from https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab, and “ViolinPlots” forked from https://github.com/bastibe/Violinplot-Matlab. <br>

Last successfully tested in MATLAB 2020b on July 19, 2023. It takes about 10-20 minutes to run the above block of code on my 2018 MacBook Pro (2.7 GHz Quad-Core Intel Core i7).  Please note this repository requires MATLAB’s Text Analytics Toolbox.  <br>

*Reference*: <br>
Koffman, B, Osman, M., Criscitiello, A., Guest, S. Women’s collaboration helps close gender gap in ice core science, *in revision*.
