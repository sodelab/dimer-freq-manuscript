
mkdir -p opt

# opt monomer
../tools/optimize/optimize.x common/monomer.xyz > opt/monomer.traj.xyz
tail -5 opt/monomer.traj.xyz > opt/monomer.opt.xyz

# opt dimer
../tools/optimize/optimize.x common/dimer.xyz > opt/dimer.traj.xyz
tail -8 opt/dimer.traj.xyz > opt/dimer.opt.xyz


