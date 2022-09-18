# VTrails Toolkit

A *vectorial* representation of the vascular network, which embodies quantitative features such as location, direction, scale, and bifurcations, has many potential cardio- and neuro-vascular applications.

VTrails stands as an end-to-end approach to extract geodesic vascular minimum spanning trees from angiographic data by solving a connectivity-optimised anisotropic level-set over a voxel-wise tensor field representing the orientation of the underlying vasculature.

VTrails has been presented at the biennial conference *Information Processing in Medical Imaging* (IPMI 2017) [1] [Full-Text](https://arxiv.org/abs/1806.03111), and further published on the *IEEE Transactions on Medical Imaging* journal [2] [Full-Text](https://ieeexplore.ieee.org/document/8421255/).

<iframe width="560" height="315" src="https://www.youtube.com/embed/ZkX3l_FhiLE" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>

## Open-Source Implementation
  
VTrails is available as Matlab toolkit.

[VTrails on YouTube](https://www.youtube.com/channel/UCC24bCFUO9uhUBLNQk1zjJw)

### Requirements:
  - Full Matlab R2016b or later
  - Correctly Configured Matlab Compiler (*Coder* and *mex* with most recent version gcc, clang, ...) 
  - Nifti Image I/O library (included) - original [source](https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
  - Boost_1_49_0 libraries (included) - original [source](https://www.boost.org/users/history/version_1_49_0.html)
  - linterp libraries (included) - original [source](http://rncarpio.github.io/linterp/)

### Installation guide:
Install VTrails Toolkit from Matlab command window by typing:
  
  `>> INSTALL_VTrails;`

### Processing Pipeline:
After successfully installing the toolkit, process the angiographic image as follows:

 > Filter with SLoGS and synthesise the Riemannian Vesselness

 > Infer the Geodesic Vascular Graph with the Connectivity-optimised Anisotropic Level-Set

 > Extract the vascular Minimum Spanning Tree(s)

   -- For more details, see `/main/PIPELINE.m` --

### Available Demos:
Two demos are included in the release:

 `>> /test_imgs/VTrailsDEMO_01;`
 
 `>> /test_imgs/VTrailsDEMO_02;`

## Developing *new* Features - Coming Soon! -

 - Topologically aware evaluation framework: quantitative accuracy assessment;
 - Alignment and registration of vascular topologies [3] [Full-Text](https://arxiv.org/abs/1809.05499);
 - Vessels morphology, lumen segmentation and clinical indices extraction;
 - Much more to come! ... Stay Tuned :)

## References
[1]
```
@InProceedings{Moriconi2017VTrails,
 title = {VTrails: Inferring Vessels with Geodesic Connectivity Trees},
 author = {Stefano Moriconi and 
           Maria A. Zuluaga and 
           Hans Rolf J{\"a}ger and 
           Parashkev Nachev and 
           S{\'e}bastien Ourselin and
           M. Jorge Cardoso},
  booktitle = {IPMI},
  year = {2017}}
```

[2]
```
@Article{Moriconi2018Inference,
 title = {Inference of Cerebrovascular Topology with Geodesic Minimum Spanning Trees},
 author = {Stefano Moriconi and 
           Maria A. Zuluaga and 
           Hans Rolf J{\"a}ger and 
           Parashkev Nachev and 
           S{\'e}bastien Ourselin and
           M. Jorge Cardoso},
 journal  = {IEEE Transactions on Medical Imaging},
 year = {2018}}
```

[3]
```
@InProceedings{Moriconi2018Elastic,
 title = {Elastic Registration of Geodesic Vascular Graphs},
 author = {Stefano Moriconi and 
           Maria A. Zuluaga and 
           Hans Rolf J{\"a}ger and 
           Parashkev Nachev and 
           S{\'e}bastien Ourselin and
           M. Jorge Cardoso},
 booktitle = {MICCAI},
 year = {2018}}
```

## Contribute to the project!

Did you like what you found here?
Then please follow the project, give us a star and help us improving the open-source code with valuable feedbacks, contributions and comments!
We will come back to you as soon as possible.

Have we made your life somehow easier? Please consider donating to support the project!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=66CWS8MB7PSGG)

We do really appreciate your generosity. Thanks for making the difference!

### Authors and Contributors
This project is currently maintained by [Stefano Moriconi](https://stefanomoriconi.github.io/mypage/) (@stefanomoriconi).
