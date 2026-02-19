# 2026-UNSW-Astronomy
Kinematic Morphology-Density Relation in MAGPI survey Clusters project. Relevant Python Scripts.

## ngist pipeline:

masterConfig.py outputs all the required .yaml files to bulk run the ngist pipeline on an entire cluster. As well as a .sh file to run all these iterations of the pipeline

- katana to run ngist pipeline 

plotKinematicMaps.py outputs a map of V (Line of Sight Stellar Velocity) across all voronoi bins, rerun any files which failed in katana and rerun the script in a terminal window rather than a batch job.
