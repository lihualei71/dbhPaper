#!/bin/bash

python3 ~/utils/pyparams.py -f "dBH_lm_params.txt" \
    'list(range(200))' \
    '[2020]' \
    '[1500, 3000]' \
    '[1000]' \
    '[0.03]' \
    '[1]' \
    '["two", "right"]' \
    '[5]' \
    '["TRUE"]' \
    '["TRUE"]'
    
