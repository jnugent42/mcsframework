# .bashrc

alias maus="source /home/ppe/j/jnugent/private/PhD_Year_1/software/MAUS-v2.9.1/env.sh"

maus 
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/lib64" 
export LD_LIBRARY_PATH="${MAUS_THIRD_PARTY}/third_party/install/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="${MAUS_THIRD_PARTY}/third_party/build/gcc-4.9.3/obj/x86_64-unknown-linux-gnu/libstdc++-v3/src/.libs:$LD_LIBRARY_PATH" 
export LD_LIBRARY_PATH="/home/ppe/j/jnugent/workarea/Unfolding/RooUnfold-1.1.1:$LD_LIBRARY_PATH"

export PYTHONPATH="/home/ppe/j/jnugent/workarea/Unfolding/Cobb_results:$PYTHONPATH"
export ROOTSYS=$MAUS_THIRD_PARTY/third_party/build/root
