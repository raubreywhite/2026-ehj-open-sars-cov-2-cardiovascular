# SSI source files

These three files are byte-for-byte copies of files inside two zip archives from Statens Serum Institut (SSI), downloaded on 25 September 2026. They are © Copyright Statens Serum Institut. The MIT licence of this repository does not cover them. `Run.R` reads all three. It stops if an md5 differs from the value below.

| File | Zip archive | Path in the zip | File md5 |
|---|---|---|---|
| `24_reinfektioner_daglig_region.csv` | https://files.ssi.dk/covid19/overvagning/dashboard/overvaagningsdata-dashboard-covid19-04012023-soc6 (zip md5 b46c633472865b8a58c93134ca342b6b) | `Regionalt_DB/24_reinfektioner_daglig_region.csv` | c3b069badfdaaf40b0d4f903d9abed26 |
| `18_fnkt_alder_uge_testede_positive_nyindlagte.csv` | the same zip as the `24_` file | `Regionalt_DB/18_fnkt_alder_uge_testede_positive_nyindlagte.csv` | 58245bd78ff14fda0831f04009329575 |
| `2026-09-23_dk_wastewater_data.csv` | https://en.ssi.dk/-/media/arkiv/uk/surveillance-and-preparedness/surveillance-in-dk/wastewater/wastewaterdata/data-wastewater-sarscov2-week38-2026-dla5.zip (zip md5 6942da2b940a7a95c84b0187ec5c8f86) | `SARS-CoV-2 koncentration/2026-09-23_dk_wastewater_data.csv` | f1ef3aff3ec7a4987c1f1a31c78fcc50 |

What `Run.R` takes from each file:

- `24_reinfektioner_daglig_region.csv` gives confirmed first infections per day and region. `Run.R` sums them per month over the 5 regions. The monthly counts from 11 March 2022 MUST equal `CFG$late_n`, which dates the recorded infections after that day. The counts from March 2020 to December 2022 are the SSI line in `figures/cohort_infection_timing.png`.
- `18_fnkt_alder_uge_testede_positive_nyindlagte.csv` gives persons tested and positive per ISO week and functional age group, as counts and per 100,000 residents. Its positives include reinfections. The confirmed cases per resident, ISO weeks 2020-W10 to 2022-W52, MUST equal `CFG$inf_age$rr` within 1e-6.
- `2026-09-23_dk_wastewater_data.csv` gives the national weekly SARS-CoV-2 RNA concentration in wastewater. Its `rna_mean_faeces` for 2022-W10 to 2022-W52 MUST equal `CFG$ww_w`, which dates the unrecorded infections after 10 March 2022.
