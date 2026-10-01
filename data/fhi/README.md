# NIPH/FHI source file

`Total forekomst_ Utvalgte diagnoser, antall pasienter.xlsx` is a byte-for-byte copy of a table from the cardiovascular disease registry statistics of the Norwegian Institute of Public Health (Folkehelseinstituttet, NIPH/FHI), downloaded on 1 October 2026. It is published under the Norwegian Licence for Open Government Data (NLOD). `Norway.R` reads it. It stops if the md5 differs from the value below.

| File | Source | File md5 |
|---|---|---|
| `Total forekomst_ Utvalgte diagnoser, antall pasienter.xlsx` | https://statistikk.fhi.no/hkr/SEZoMIr4ryq7rCm076P4WnYbXmhJK3X8?Argang=2016,2017,2018,2019,2020,2021,2022,2023,2024,2025&Kjonn=Alle&Alder_FireDelt=Alle,0-49,50-69,70-89,90%2B&Kilde=Alle&Gruppe=Pasienter,Sirkulasjonssystemet,Hypertensjon,IskemiskHjertesykdom,AnginaPectoris,Hjerteinfarkt,AtrieflimmerAtrieflutter,Hjerneslag,Hjertesvikt,TIA,Brystsmerter&MEASURE_TYPE=Antall | 6750b08e2fd247555b2f35d061a14475 |

The table gives, per year from 2016 to 2025, the number of patients with at least one contact for each diagnosis group ("total forekomst"), both sexes, all sources, for all ages and for the age groups 0-49, 50-69, 70-89 and 90 or more. `Norway.R` uses the 4 age groups for all CVD (I00-I99), ischemic heart disease, atrial fibrillation and flutter, stroke and heart failure. It checks that the age groups sum to the all-ages row.
