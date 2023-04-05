# kpf_cap1h
Code to compute the throughput of KPF from spectrophotometric observations

Some notes:

- calc_throughput_tools.py contains useful functions for computing the throughput

- examples/ contains example scripts that uses the aforementioned functions to compute and plot the throughput for KPF main, the exposure meter, and the H&K spectrometer (modes = "kpf", "em", and "hk" respectively)

- inputs/ contains the data used in the computation - the telluric file needs to be unzipped and the KPF data needs to be downloaded. The file lists in the stellar folders list the KPF data associated with spectrophotometric observations for that star and includes data where the telescope pointing was offset. The example scripts are setup with the file name of the observation that showed the best coupling
