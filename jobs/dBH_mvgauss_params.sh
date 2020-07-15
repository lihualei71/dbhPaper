#!/bin/bash

python3 ~/utils/pyparams.py -f "dBH_mvgauss_params.txt" \
    'list(range(100))' \
    '[1000]' \
    '[0.01]' \
    '[1]' \
    '["two", "right"]' \
    '[10]' \
    '["TRUE"]'
    
