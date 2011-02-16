BatchSystems = {
   PBS = {
      launchpad = {
         submitHeader= [[
               #PBS -N tm-$(JOBNAME)
               #PBS -l walltime="$(TIME)"
               #PBS -l nodes=1:GPU:ppn=8
               #PBS -o std.out
               #PBS -e std.err
               #PBS -V
               #PBS -M rge21@nmr.mgh.harvard.edu
               #PBS -m a
               cd $PBS_O_WORKDIR
         ]],
         maxCoresPerNode = 8,
      },
   },
}
