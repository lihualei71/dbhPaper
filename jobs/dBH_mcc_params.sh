#!/bin/bash

python3 ~/utils/pyparams.py -f "dBH_mcc_params.txt" \
    'list(range(100))' \
    '[100]' \
    '[3, 30]' \
    '[0.3]' \
    '[1]' \
    '["two", "right"]' \
    '[10]' \
    '["TRUE"]' \
    '["TRUE"]'
    
