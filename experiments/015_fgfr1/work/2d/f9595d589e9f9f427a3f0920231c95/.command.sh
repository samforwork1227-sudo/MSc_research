#!/bin/bash -euo pipefail
start_ts=$(date +%s)
    echo "[RUNCHIMERAX] Started at $(date)"

    export HOME=$PWD
    export XDG_CONFIG_HOME=$PWD/.config
    export XDG_DATA_HOME=$PWD/.local/share
    export XDG_CACHE_HOME=$PWD/.cache
    mkdir -p "$HOME" "$XDG_CONFIG_HOME" "$XDG_DATA_HOME" "$XDG_CACHE_HOME"

    echo "PWD=$PWD"
    echo "ciffile=fgfr1_20lipids_model_0.cif"
    ls -l
    ls -l "fgfr1_20lipids_model_0.cif"

    cat > run.cxc <<EOF
close session
open fgfr1_20lipids_model_0.cif
cartoon
color /A:1-376 blue
color /A:377-397 red
color /A:398-822 green

show #1/B-U
style #1/B-U stick
color #1/B-U gold

measure inertia #1/A:377-397@CA

shape cylinder fromPoint 8.6225,-11.792,-4.3179 axis -0.422,0.112,-0.900 height 40 radius 0.3 color red coordinateSystem #1
shape cylinder fromPoint 8.6225,-11.792,-4.3179 axis 0.422,-0.112,0.900 height 40 radius 0.3 color red coordinateSystem #1

# ICD
torsion #1:A:408@C #1:A:409@N #1:A:409@CA #1:A:409@C 150 move large
torsion #1:A:419@C #1:A:420@N #1:A:420@CA #1:A:420@C -8 move small
torsion #1:A:399@C #1:A:400@N #1:A:400@CA #1:A:400@C 10 move large

# ECD hinges
torsion #1:A:253@C #1:A:254@N #1:A:254@CA #1:A:254@C -10 move small
torsion #1:A:376@C #1:A:377@N #1:A:377@CA #1:A:377@C 30 move small
torsion #1:A:154@C #1:A:155@N #1:A:155@CA #1:A:155@C 15 move small
torsion #1:A:376@C #1:A:377@N #1:A:377@CA #1:A:377@C 10 move small
torsion #1:A:376@C #1:A:377@N #1:A:377@CA #1:A:377@C 55 move small
torsion #1:A:250@C #1:A:251@N #1:A:251@CA #1:A:251@C 12 move small
torsion #1:A:255@C #1:A:256@N #1:A:256@CA #1:A:256@C 25 move small
torsion #1:A:360@C #1:A:361@N #1:A:361@CA #1:A:361@C 20 move small
torsion #1:A:370@C #1:A:371@N #1:A:371@CA #1:A:371@C 30 move small

save "$PWD/protein_preMD.pdb" format pdb models #1
exit
EOF

    echo "===== run.cxc ====="
    cat run.cxc

    chimerax --nogui --script "$PWD/run.cxc"
    echo "ChimeraX exit code: $?"

    ls -l
    test -s protein_preMD.pdb

    end_ts=$(date +%s)
    elapsed=$(( end_ts - start_ts ))
    printf "[RUNCHIMERAX] Finished in %d min %d s\n"       $((elapsed/60)) $((elapsed%60))
