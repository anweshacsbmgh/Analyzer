logout
logoff
quit
logout
MATLAB -nojit
bsub -Is -q interactive bash
bsub -Is -q interactive -R "rusage[mem=16000]" bash
module avail
module load math/matlab/2015a
which matlab
matlab nosplash -nojit
bsub -q short -w 12:00 matlab autoreadAllChangi.m
bsub -q short -W 12:00 matlab autoreadAllChangi.m
