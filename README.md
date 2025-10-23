# Wind-Stress-and-Wind-Stress-Curl-Calculation

This repository provides Python implementations for calculating **wind stress** and **wind stress curl**, offering an accessible alternative to existing MATLAB and FORTRAN scripts.  
Building upon **Ramkrushn Patel’s MATLAB code**, this version refines the **drag coefficient** for high wind speeds and corrects several minor issues to improve accuracy and usability.

The drag coefficient is calculated as follows: <br>

|Cd| Wind Speed |
|---|---|
| 0.00218 | (U <= 1) |
| (0.62 + 1.56 / U) * 0.001|(1 < U <= 3)|
| 0.00114|(3 < U <= 10)|
| (0.49 + 0.065 * U) * 0.001|(10 < U <= 19)|
| (1.364 + 0.0234 * U - 0.0002 * U ** 2) * 0.001 |(U > 19)|


<!-- $$ 
C_d=
\begin{cases}
0.00218 & ,U <= 1 \\
(0.62 + 1.56 / U) \times 0.001 & ,1 < U <= 3 \\
0.00114 & ,3 < U <= 10 \\
(0.49 + 0.065 \times U) \times 0.001 & ,10 < U <= 19 \\
(1.364 + 0.0234 \times U - 0.0002 \times U^2) \times 0.001 & ,U > 19
\end{cases}
$$ -->

![C_d formula](https://latex.codecogs.com/svg.image?C_d=%20%5Cbegin%7Bcases%7D%200.00218%20%26%20,U%20%5Cleq%201%20%5C%5C%20(0.62%20%2B%201.56%20%2F%20U)%20%5Ctimes%200.001%20%26%20,1%20%3C%20U%20%5Cleq%203%20%5C%5C%200.00114%20%26%20,3%20%3C%20U%20%5Cleq%2010%20%5C%5C%20(0.49%20%2B%200.065%20%5Ctimes%20U)%20%5Ctimes%200.001%20%26%20,10%20%3C%20U%20%5Cleq%2019%20%5C%5C%20(1.364%20%2B%200.0234%20%5Ctimes%20U%20-%200.0002%20%5Ctimes%20U%5E2)%20%5Ctimes%200.001%20%26%20,U%20%3E%2019%20%5Cend%7Bcases%7D)

where the $U$ is the value of wind speed in $m/s$.

---

## 🌀 Usage

```python
from windstress import ra_windstr
from windstress_curl import ra_windstrcurl

# Calculate wind stress components
taux_cal, tauy_cal = ra_windstr(u, v)

# Calculate wind stress curl
WSC_cal = ra_windstrcurl(lat, lon, u, v)
```
---

## 📊 Validation
### Compare our results with the initial wind stress and wind stress curl data
We randomly choose a snapshot (2025-10-20T 12:00) from [Global Ocean Hourly Sea Surface Wind and Stress from Scatterometer and Model](https://data.marine.copernicus.eu/product/WIND_GLO_PHY_L4_NRT_012_004/description) data to validate the our scripts' result.
We use the initial wind to compute the wind stress and wind stress curl. Our result is very simillar to the raw data. 

### **Figure 1.** Comparison of Wind Stress
![Results: Wind Stress](output.png)

### **Figure 2.** Comparsion of Wind Stress Curl
![Results: Wind Stress Curl](output2.png)

---

## 📁 Repository Structure
```shell
Wind-Stress-and-Wind-Stress-Curl-Calculation/
├── windstress.py                # Compute wind stress from wind components
├── windstress_curl.py           # Compute wind stress curl
├── test_compare.ipynb           # Validation and visualization notebook
├── cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H_1761115052002.nc  # CMEMS wind data
├── output.png                   # Wind stress comparison figure
├── output2.png                  # Wind stress curl comparison figure
├── LICENSE
└── README.md
```

