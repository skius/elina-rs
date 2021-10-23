#!/bin/bash

bindgen wrapper.h --blocklist-item "FP_NAN"\
                  --blocklist-item "FP_INFINITE"\
                  --blocklist-item "FP_ZERO"\
                  --blocklist-item "FP_SUBNORMAL"\
                  --blocklist-item "FP_NORMAL"\
                  -o src/bindings.rs