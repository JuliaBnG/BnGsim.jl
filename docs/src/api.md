# API Reference

```@autodocs
Modules = [BnGsim, BnGsim.Domain, BnGsim.ExitCode, BnGsim.FisherWrightAdapter]
```

## Legacy building blocks (0.1–0.3)

`BnGsim.Legacy` holds the building blocks of BnGsim 0.1–0.3. They draw
from an explicit `rng` or the global generator, not from the keyed
streams of the breeding loop, so they are outside the reproducibility
invariants and the 1.0 promise. Load them with `using BnGsim.Legacy`.

```@autodocs
Modules = [BnGsim.Legacy]
```
