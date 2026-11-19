# Flight and weather data

The flight records and several weather files were obtained from the RocketPy
repository. This project uses those input files; it does not run RocketPy.
JSBSim performs the flight simulations, and Python scripts process the records.
Original input files are retained alongside the derived profiles and results.

## Flight records and vehicle inputs

| Files | Source and use |
| --- | --- |
| `val/data/2022-06-24-serial-5115-flight-0001-TeleMetrum.csv` and `2022-06-24-serial-6583-flight-0003-TeleMega.csv` | Western Engineering Rocketry Team, Prometheus flight on 24 June 2022. The RocketPy repository identifies the vehicle under `data/rockets/prometheus`. These are the GPS, pressure and barometric-height records used in the comparison. |
| `val/data/Cesaroni_7579M1520-P.eng` | Prometheus motor thrust curve supplied with the vehicle inputs. |
| `val/data/cavour/altimeter_cavour.csv` | PoliTo Rocket Team, Cavour flight at Ponte de Sor on 13 October 2023. The RocketPy repository identifies the vehicle under `data/rockets/polito`. This is a barometric-height record. |
| `val/data/cavour/Cesaroni_3618L995-P.eng` | Cavour motor thrust curve supplied with the vehicle inputs. |
| `val/data/cavour/drag_coefficient_power_off.csv` and `drag_coefficient_power_on.csv` | Supplied Cavour aerodynamic tables. The model uses the power-off table during both boost and coast. Both source tables are retained. |
| `val/data/cavour/euroc_2023_all_windows.nc` | Included ERA5 weather data used to interpolate the Cavour launch-time profile between 13:00 and 14:00 UTC. |

For Cavour, the simulated geometric apogee is converted to pressure height
using a standard dry ICAO atmosphere, zeroed at the pad, before comparison
with the barometric record. The original instrument's conversion is not
documented, so this is the assumed comparison convention. The resulting
apogee discrepancy is -2.03%.

## Prometheus regional weather

`val/data/prometheus_2022_igra_excerpt.txt` contains the two original NOAA IGRA
records for station USM00072364 at 12:00 UTC on 24 June and 00:00 UTC on
25 June 2022. The station is about 120 km from the launch pad.
The full station archive was obtained from
<https://www.ncei.noaa.gov/pub/data/igra/data/data-por/USM00072364-data.txt.zip>.
Its hash and the relationship between that archive and the included excerpt
are recorded in `provenance/data_sources.json`.

`provenance/igra2-data-format.txt` and `igra2-station-list.txt` retain NOAA's
format and station information. `val/prepare_prometheus_weather.py` produces
the two `prometheus_sounding_*.csv` files and `prometheus_inputs_2022.json`.
Valid temperature/pressure levels and valid wind levels are selected separately.
Missing wind therefore does not remove an otherwise valid temperature record.

The input metadata records a 19.20 K density fit over 1,401–6,101 m geometric
altitude above sea level. The flight builder refits the shorter 1,401–5,701 m
column and uses 19.25 K. These are two fit ranges, not two weather datasets.

## Regional, flight-informed and proxy comparisons

The regional-weather Prometheus prediction has a +3.08% geometric-apogee error
against the fitted GPS reference. The flight-informed comparison uses pressure
and GPS from that same flight to infer a 16.74 K column-temperature offset.
With the same 2022 wind inputs, its geometric-apogee error is +2.83%.
This second comparison is calibration using the flight, not independent weather
validation. Pressure-selection sensitivity is documented separately in
`rq2_revision2/prometheus_temperature/pressure_selection_diagnostic.json`.

The 2023 proxy-weather case is retained to explain the earlier +1.29% result.
Its files include `val/prometheus_inputs.py`, `val/prometheus_inputs.json` and
`val/data/spaceport_america_pressure_levels_2023_hourly.nc`.
Although that JSON contains the proxy profile, the current weather preparer
also reads it for the vehicle's base inputs, so it must remain in place.
The current flight builder defaults to `prometheus_inputs_2022.json`.

## ERA5 files and source records

The three original netCDF files in `atmos/` support the column-temperature and
persistence statistics. They describe the selected sites and periods; they do
not directly identify an uncertainty distribution for every launch location.
`provenance/data_sources.json` records their original URLs, SHA-256 hashes and
RocketPy commit `9bd6ad3af8f97bafa3201d4e70eba2877e66e040`.
It also identifies the 2023 Spaceport proxy file at that commit.

The later temperature correction is recorded in
`rq2_revision2/prometheus_temperature/supplement_acceptance.json` and
`rq2_revision2/prometheus_temperature/source_integration_proof.json`; these supplement the earlier
`results_revised/validation/validation_run_manifest.json` for the changed
temperature diagnostics while retaining the original two 200-case studies
(Prometheus seed 7 and Cavour seed 11).
Historical records retain their original bytes and execution paths.
