#!/bin/bash -euo pipefail
start_ts=$(date +%s)
    echo "[RUNCHIMERAX] Started at $(date)"

    export HOME="$PWD"
    export XDG_CONFIG_HOME="$PWD/.config"
    export XDG_DATA_HOME="$PWD/.local/share"
    export XDG_CACHE_HOME="$PWD/.cache"
    mkdir -p "$HOME" "$XDG_CONFIG_HOME" "$XDG_DATA_HOME" "$XDG_CACHE_HOME"

    cat > run.cxc <<EOF
    # subdomains & linkers
    #1:A:25-125      # sub1
    #1:A:126-153     # linker12
    #1:A:154-247     # sub2
    #1:A:248-255     # linker23
    #1:A:256-358     # sub3
    #1:A:481-670     # ICD
    #1:A:378-398     # TM

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

    shape cylinder fromPoint 7.4935,11.63,9.6412 axis -0.582,0.739,-0.338 height 40 radius 0.3 color red coordinateSystem #1
    shape cylinder fromPoint 7.4935,11.63,9.6412 axis 0.582,-0.739,0.338 height 40 radius 0.3 color red coordinateSystem #1

    # ICD
    torsion #1:A:408@C #1:A:409@N #1:A:409@CA #1:A:409@C 90 move large
    torsion #1:A:399@C #1:A:400@N #1:A:400@CA #1:A:400@C 150 move large
    torsion #1:A:470@C #1:A:471@N #1:A:471@CA #1:A:471@C 10 move small
    torsion #1:A:419@C #1:A:420@N #1:A:420@CA #1:A:420@C -8 move small
    torsion #1:A:420@C #1:A:421@N #1:A:421@CA #1:A:421@C -10 move small
    torsion #1:A:418@C #1:A:419@N #1:A:419@CA #1:A:419@C 10 move small

    getcrd #1:A:401@CA
    turn 0,1,0 -60 atoms #1:A:401-412 center 18.061,-0.777,8.716 coordinateSystem #1

    # ECD
    turn 0.582,-0.739,0.338 8 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:361-377
    turn 0.582,-0.739,0.338 10 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:256-360
    turn 0.582,-0.739,0.338 10 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:256-360
    turn -0.582,0.739,-0.338 5 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:256-360

    turn 0.582,-0.739,0.338 12 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:154-255
    turn 0.582,-0.739,0.338 8 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:154-255
    turn 0.582,-0.739,0.338 10 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:25-153

    torsion #1/A:250@C #1/A:251@N #1/A:251@CA #1/A:251@C 10 move small
    torsion #1/A:253@C #1/A:254@N #1/A:254@CA #1/A:254@C -10 move small

    torsion #1/A:376@C #1/A:377@N #1/A:377@CA #1/A:377@C 30 move small

    torsion #1/A:152@C #1/A:153@N #1/A:153@CA #1/A:153@C 20 move small
    torsion #1/A:154@C #1/A:155@N #1/A:155@CA #1/A:155@C 15 move small

    save "protein_preMD.pdb" format pdb models #1
    exit
    EOF

    echo "Workdir contents before ChimeraX:"
    ls -l

    chimerax --nogui --script run.cxc
    echo "ChimeraX exit code: $?"

    ls -l
    test -s protein_preMD.pdb

    end_ts=$(date +%s)
    elapsed=$(( end_ts - start_ts ))
    printf "[RUNCHIMERAX] Finished in %d min %d s
"       $((elapsed/60)) $((elapsed%60))
