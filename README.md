# FSTSP-DR-LRP

C++17 implementation of the two-phase matheuristic for the flying sidekick traveling salesman problem with drone resupply (FSTSP-DR-LRP). The solver returns a JSON timetable with the truck route, delivery-drone sorties and resupply-drone sorties.

## Historical Solomon-derived inputs

`data/solomon/` contains the 144 archived input files associated with the main-text experimental design. The file names encode the historical distribution label and the three input-generation factors. All four label families in this archive share the C101 coordinates. Use the labels as archived identifiers. The `DELAY_TIME` column stores relative release coefficients, with blank entries for initially available parcels. The absolute release scale used for each historical run is a separate input.

Convert an input with Python 3 using `python prepare_solomon.py SOURCE.txt INSTANCE.txt --release-scale SCALE --rendezvous-count COUNT`. Supply the scale and number of candidate stops for your run. The converter uses the closest customers to the depot, breaking ties by customer ID. Speeds, capacity and endurance in the generated instance are physical experiment settings.

## Build and run

Use a C++17 compiler, CMake and the Gurobi 12 C++ SDK. Set `GUROBI_HOME` to the SDK directory and make its shared libraries available to the executable.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Copy `algorithm.params` to a local parameter file and enter values for its fields. Run `fstsp_dr_lrp INSTANCE.txt --config LOCAL_PARAMS`. Seed and time limit can also be supplied with `--seed` and `--seconds`. Command-line values take precedence. The JSON output includes the effective search parameters for that run.
