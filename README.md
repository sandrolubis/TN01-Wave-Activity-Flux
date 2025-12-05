# Takaya–Nakamura Wave-Activity Flux (NCL)

This repository provides an NCL implementation of the quasi-geostrophic wave-activity flux formulated by Takaya & Nakamura (1997, 2001).  
This script computes horizontal wave-activity flux components (Fx, Fy) and the QG streamfunction anomaly (ψ′) at a given pressure level.

By **Sandro W. Lubis** (PNNL). The code has been used in several studies of large-scale teleconnections and extreme events, including:

Lubis, S. W., Chen, Z., Lu, J., Hagos, S., Chang, C.-C., & Leung, L. R. (2024).
Enhanced Pacific Northwest heat extremes and wildfire risks induced by the boreal summer intraseasonal oscillation.
npj Climate and Atmospheric Science. https://doi.org/10.1038/s41612-024-00766-3

---

## Features

- Implements Takaya–Nakamura (2001) Eq. (38) for stationary/migratory Rossby waves.
- Uses:
  - Daily geopotential height anomalies (e.g., ERA5) at one pressure level (default: 200 hPa).
  - Monthly/seasonal climatological U, V, and Z as basic state.
- Outputs:
  - `Fx` – zonal component of wave-activity flux (m² s⁻²)
  - `Fy` – meridional component of wave-activity flux (m² s⁻²)
  - `psidev` – QG streamfunction anomaly (m² s⁻¹)
- Loops over multiple years and writes yearly NetCDF files.

---

## Requirements

- [NCL (NCAR Command Language)](https://www.ncl.ucar.edu/)
- NCAR libraries:
  - `gsn_code.ncl`
  - `gsn_csm.ncl`
  - `contributed.ncl`
- NetCDF I/O support

Input data:

- Daily geopotential height (or geopotential) anomalies: `z3.YYYY.nc`
- Monthly/seasonal climatologies at same level:
  - Zonal wind (U)
  - Meridional wind (V)
  - Geopotential height (Z)

---

## Directory Structure

Example structure:

```text
.
├── tn2001_waf.ncl           # Main script
└── input
    └── climo
        ├── u200.seas.nc     # U climatology (var131)
        ├── v200.seas.nc     # V climatology (var132)
        └── z200.seas.nc     # Z climatology (var129)
