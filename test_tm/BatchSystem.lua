BatchSystems = {
   PBS = {
      launchpad = {
         submitHeader= [[
               #PBS -N tm-$(JOBNAME)
               #PBS -l walltime="$(TIME)"
               #PBS -l nodes=1:ppn=8:GPU
               #PBS -o std.out
               #PBS -e std.err
               #PBS -V
               #PBS -q GPU
               #PBS -M rge21@nmr.mgh.harvard.edu
               #PBS -m a
               cd $PBS_O_WORKDIR
               echo $HOSTNAME
         ]],
         maxCoresPerNode = 8,
      },
   },
}
