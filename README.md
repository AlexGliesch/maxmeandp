# A VLNS heuristic for the Max-Mean Dispersion Problem

Toy implementation of a VLNS heuristic for the Max-Mean Dispersion Problem (MMDP), which I made because I wanted to learn Rust. The heuristic works ok and is actually better than the VNS of [Brimberg et al. (2017)](http://dx.doi.org/10.1016/j.ins.2016.12.021), but worse than the memetic algorithm of [Lai et al. (2020)](https://doi.org/10.1016/j.eswa.2019.112856). For a description of the MMMDP, see the aforementioned papers.

A full description of the heuristic and test results are coming.

## Running the code

1. Unpack Brimberg et al. (2017)'s instances with `tar -zxvf instances-brimberg.tar.gz` (~1.2GB).
1. Generate Lai et al. (2020)'s instances running `./generate-lai-instances.sh` under `resources/lai` (~26GB).
1. Build using `cargo --build release`.
1. Run using `cargo run -- -i {instance} -t {timeLimit} -s {seed}`. For more options, see `--help`.
