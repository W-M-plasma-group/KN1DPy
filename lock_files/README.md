## Lock files

`pixi.lock` in this folder is the exact environment used to generate the reference outputs in `examples/`. Use it if you want to guarantee bit-for-bit reproduction of those results.

**To use it:**

1. Copy `pixi.lock` from this folder to the KN1DPy root directory
2. Run `pixi install` — pixi will use the lock file instead of re-solving