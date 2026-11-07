# A lightweight allocation/time probe. Run after verifying correctness.
using MMIEF
p=Config(rmax=80.0,dr=0.25,tmax=1.0,pulse_amplitude=0.0)
r,u=MMIEF.InitialData.initial_data(p)
w=MMIEF.Evolution.Workspace(r)
MMIEF.Evolution.setup_vacuum_balance!(p,w)
du=similar(u)
MMIEF.Evolution.rhs!(du,u,0.0,p,w)
println("RHS allocations (bytes): ",@allocated MMIEF.Evolution.rhs!(du,u,0.0,p,w))
println("RHS elapsed time (seconds): ",@elapsed MMIEF.Evolution.rhs!(du,u,0.0,p,w))
