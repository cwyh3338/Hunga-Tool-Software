# Hunga-Tool-Software

The Hunga Tool retrieves stratospheric volcanic aerosols and SO2 using hyperspectral solar backscattered UV (BUV) radiance and radiative transfer modeling.
Repository: [https://github.com/cwyh3338/Hunga-Tool-Software](https://github.com/cwyh3338/Hunga-Tool-Software)

---

## Features
- Supports 2-parameter retrieval (AOD, peak height) and 3-parameter retrieval (AOD, peak height, SO2 column).
- Uses VLIDORT v2.8 radiative transfer model and Mie scattering calculations.
- Input data include hyperspectral BUV radiance (e.g., TROPOMI, OMPS, OMI) and atmospheric profiles (temperature, pressure, ozone).
- Output includes retrieved AOD, plume height, SO2 VCD, fitting diagnostics, and Jacobians.

---

## Current Version of The Package
- HTHH_RETRIEVALS_V4_with_SO2_28oct2024_GSFC 
  (Package version dated October 28, 2024)

---

## 📌 Note to Readers (August 28, 2025)
This retrieval tool is available here for you to look at and run if you wish.  
However, please note that the code is still undergoing refinement as we work on the follow-on paper. Planned updates include:
- a simple treatment for tropospheric clouds,
- minor bug fixes,
- improvements to make certain inputs more user-friendly.

The current version does not yet include a treatment for tropospheric cloud interference in the extended fitting window required for 3-parameter retrievals. For this reason, we recommend using the 2-parameter retrievals (AOD, peak height) with wavelengths below ~296 nm in the present release.

You are welcome to try the current version, but you may prefer to wait for the updated release.


---

## Documentation
For detailed descriptions about structure of the package, how to run and prepare input files, and output files, please refer to the uploaded user guide file ([HTHH User's Guide] Aug28 2025.pdf)

---

## Installation

Clone this repository and move into the package directory:

```bash
git clone https://github.com/cwyh3338/Hunga-Tool-Software.git
cd Hunga-Tool-Software/HTHH_RETRIEVALS_V4_with_SO2_28oct2024_GSFC
```

---

## Reference
Spurr, R. J. D., et al. (2025): Solar Backscatter Ultraviolet (BUV) Retrievals of Mid-Stratospheric Aerosols from the 2022 Hunga Eruption, EGUsphere [preprint]. DOI: 10.5194/egusphere-2025-2938

---

## Contact
Wonei Choi – won-ei.choi@nasa.gov

