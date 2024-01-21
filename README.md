# DiaphragmTracking
Official code for the diaphragm tracking algorithm

Following the workflow depicted below:
* segmentDiaphragm.m performs diaphragm segmentation using the peak-exhale 4D-CT
* getDph3DShift.m computes the principal motion vector by rigid registration to the peak-inhale 4D-CT
* get2DDiaphragmModel.m forward projects the 3D diaphragm model to the 2D projection space
* trackDiaphragm.m tracks the diaphragm on each candidate 2D projection using a maximal gradient algorithm

For details please see Hindley, N., Keall, P., Booth, J. and Shieh, C.‐C. (2019), Real‐time direct diaphragm tracking using kV imaging on a standard linear accelerator. Med. Phys., 46: 4481-4489. doi:10.1002/mp.13738

![Proposed clinical workflow for the diaphragm tracking algorithm](https://github.sydney.edu.au/Image-X/DiaphragmTracking/blob/master/Workflow.jpg?raw=true)
