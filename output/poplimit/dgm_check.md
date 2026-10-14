# DGM check (Stage 1, step 1; protocol v1.2)
Generated 2026-10-10 10:26:53.975816

## Table 1. Marginal distributions of X1-X3 in the two trials (AC mean 50 / SD 10; BC mean 53 / SD 12; D4 truncated at 70)

| D | trial | type | mean | sd | skew | Q1 | median | Q3 | p_gt_cut | comp_ratio | mix_p | mix_m1 | mix_m2 | mix_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D1 | AC | normal | 50.000 | 10.00 | 0.000 | 43.255 | 50.000 | 56.745 | 0.067 | NA | NA | NA | NA | NA |
| D1 | BC | normal | 53.000 | 12.00 | 0.000 | 44.906 | 53.000 | 61.094 | 0.159 | NA | NA | NA | NA | NA |
| D2 | AC | mixture | 50.000 | 10.00 | 0.800 | 43.402 | 48.843 | 54.915 | 0.083 | 0.75 | 0.905 | 47.858 | 70.427 | 7.5 |
| D2 | BC | mixture | 53.000 | 12.00 | 0.800 | 45.082 | 51.611 | 58.898 | 0.135 | 0.75 | 0.905 | 50.430 | 77.513 | 9.0 |
| D3 | AC | mixture | 50.000 | 10.00 | -0.800 | 45.085 | 51.157 | 56.598 | 0.039 | 0.75 | 0.905 | 52.142 | 29.573 | 7.5 |
| D3 | BC | mixture | 53.000 | 12.00 | -0.800 | 47.102 | 54.389 | 60.918 | 0.133 | 0.75 | 0.905 | 55.570 | 28.487 | 9.0 |
| D4 | AC | normal | 50.000 | 10.00 | 0.000 | 43.255 | 50.000 | 56.745 | 0.067 | NA | NA | NA | NA | NA |
| D4 | BC | truncnorm | 53.002 | 9.95 | -0.511 | 46.424 | 53.900 | 60.700 | 0.115 | NA | NA | NA | NA | NA |

## Table 1b. Oracle-weight ESS/N (1-D and 3-D, rho = 0.35)

| D | ess_oracle_1d | ess_oracle_3d | w_q999_3d |
|---|---|---|---|
| D1 | 0.765 | 0.563 | 10.104 |
| D2 | 0.685 | 0.434 | 10.322 |
| D3 | 0.638 | 0.407 | 13.947 |
| D4 | 0.801 | 0.496 | 10.750 |

## Table 2. b_EM (shared by X1-X3) and unadjusted transport gap (marginal logOR, AC pop - BC pop); gap2 = strong EM (x2)

| shape | b_em | gap_D1 | gap_D2 | gap_D3 | gap_D4 | b_em_x2 | gap2_D1 | gap2_D2 | gap2_D3 | gap2_D4 |
|---|---|---|---|---|---|---|---|---|---|---|
| linear | 0.4139 | -0.2 | -0.182 | -0.232 | -0.199 | 0.8278 | -0.241 | -0.233 | -0.286 | -0.293 |
| quad | 0.1378 | -0.2 | -0.146 | -0.222 | -0.042 | 0.2755 | -0.304 | -0.192 | -0.360 | -0.075 |
| threshold | 0.7998 | -0.2 | -0.099 | -0.227 | -0.091 | 1.5996 | -0.297 | -0.135 | -0.343 | -0.143 |

## Table 3. Centring constants E_AC[g(X)]

| D | linear | quad | threshold |
|---|---|---|---|
| D1 | 0 | 1 | 0.0667 |
| D2 | 0 | 1 | 0.0827 |
| D3 | 0 | 1 | 0.0391 |
| D4 | 0 | 1 | 0.0667 |

Censoring rate (exponential, shared by both trials): 0.0075
