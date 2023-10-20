**EWALD** is an open-source software that implements 3 optical diffraction tomography reconstruction algorithms, all based on Fourier Diffraction Theorem (FDT): (1) DI (Direct Inversion), GP (Gerchberg-Papoulis algorithm) and GPSC (Gerchberg-Papoulis with total-variation-based object-support constraint). It is dedicated for investigation of biological and technical samples, like cells, bacterias, optical fibers etc.

EWALD is a flexible software that allows:
- reconstruction of data in **transmission** and **reflection** mode, or both
- reconstruction of data captured with **single** or **multiple wavelengths**
- reconstruction of data captured in **"limited angle" configuration** (angular scanning of the laser beam with stationary sample and camera)
- reconstruction of data captured with stationary laser beam and a camera, **with rotating sample**

Example of tomographic reconstructions can be found [here](https://biophase.pl/ewald/).

### Installation

Currently, only the Matlab version of the software is available. To get the code running:
1. Download the repository
2. Download an example measurement dataset: [Sinogram_03_hacat_03.mat](https://wutwaw-my.sharepoint.com/:u:/g/personal/wojciech_krauze_pw_edu_pl/Ed-8DABhSV1Elr5gUhCPy6IBOy4L5YOw-q1TZgDYVA5IBw?e=BtCY06).
3. Run the RecGTVIC.m file in Matlab.

If you want to use the GPSC algorithm with advanced total-variation minimization you need to make these 2 additional steps before running RecGTVIC.m:
1. Download ASTRA Tomography Toolbox [astra-toolbox/astra-toolbox: ASTRA Tomography Toolbox (github.com)](https://github.com/astra-toolbox/astra-toolbox) and add it to Matlab path.
2. Download linear-operator toolbox Spot [mpf/spot: A linear-operator toolbox for Matlab (github.com)](https://github.com/mpf/spot) and add it to Matlab path.

Note, that you need to have CUDA on your computer to run GPSC!

### Licensing

The code is shared under GPLv3 license. If you use this code, please cite one of the following papers, depending on your application:

- if you use DI or GP reconstruction methods in transmission, please cite:

[1] W. Krauze, P. Makowski, M. Kujawińska, and A. Kuś, “Generalized total variation iterative constraint strategy in limited angle optical diffraction tomography,” Opt. Express 24(5), 4924–4936 (2016).

- if you use GPSC in transmission, please cite:

[2] W. Krauze "Optical diffraction tomography with finite object support for the minimization of missing cone artifacts," Biomed. Opt. Express 11(4), 1919-1926 (2020).

- if you use multiwavelength modality, please cite:

[3] P. Ossowski et al. "Near-infrared, wavelength, and illumination scanning holographic tomography," Biomed. Opt. Express 13(11), 5971-5988 (2022).

- if you use the reflection mode, please cite:

[4] W. Krauze, P. Ossowski, M. Nowakowski, M. Szkulmowski, M. Kujawińska "Enhanced QPI functionality by combining OCT and ODT methods," Proc. SPIE 11653, 19-24 (2021).

**Main contributors:**
- Piotr L. Makowski (2014-2018)
- Paweł Ossowski (2020-2022)
- Wojciech Krauze (2016-now)

**Other contributors:** 
- Michał Ziemczonok 
- Piotr Stępień.

### Dependencies
ASTRA Tomography Toolbox [astra-toolbox/astra-toolbox: ASTRA Tomography Toolbox (github.com)](https://github.com/astra-toolbox/astra-toolbox)
pyl1 [3cHeLoN/pyl1 (github.com)](https://github.com/3cHeLoN/pyl1)
