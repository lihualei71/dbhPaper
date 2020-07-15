#!/bin/bash

python3 ~/utils/pyparams.py -f "dBH_mvt_params.txt" \
    'list(range(100))' \
    '[1000]' \
    '[50]' \
    '[0.01]' \
    '[1]' \
    '["two", "right"]' \
    '[10]' \
    '["TRUE"]'
    
python3 ~/utils/pyparams.py -f "dBH_mvt_params.txt" -a "True"\
    'list(range(100))' \
    '[100]' \
    '[5]' \
    '[0.1]' \
    '[1]' \
    '["two", "right"]' \
    '[10]' \
    '["TRUE"]'
    
