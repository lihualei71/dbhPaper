#!/bin/bash

python3 ~/utils/pyparams.py -f "dBH_mvt_timing_params.txt" \
    'list(range(100))' \
    '[10, 30]' \
    '[50]' \
    '[0.8]' \
    '[1]' \
    '[1]' \
    '["TRUE"]'
    
