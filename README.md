# HIV-1 In-Host Dynamics Simulator

A study-grade mathematical model of within-host HIV-1 dynamics in the
**stable Nowak–May form**, extended with the **latent reservoir** and a
**two-strain (R5/X4) coreceptor-escape** term.

Why two strains? It is the reason single-target CCR5 knock-out therapies keep
failing clinically: block CCR5 and the virus switches to CXCR4. The model lets
you test which interventions act on both strains.

**Interventions modelled:** ART (RTI + PI), SHOCK (latency-reversing pulses),
CCR5 gene-therapy protection, CAR-T (armed T-cells), bNAb (neutralizing window).

**State:** `x` uninfected CD4, `y5`/`y4` R5/X4-infected, `l` latent reservoir,
`v5`/`v4` free virus, `D` cumulative immune damage.

## Files
| file | purpose |
|------|---------|
| `model.py` | the ODE model + interventions |
| `miracle_hunt.py`, `champion_sweep.py`, `frontier.py`, `panlandscape.py` | parameter searches |
| `robustness50.py`, `stoch_demo.py`, `stoch_tail.py` | robustness + stochastic runs |
| `plot_*.py`, `clones.py` | figures |
| `MIRACLE_CURE.md`, `RESULTS_GO.md` | write-ups of the findings |

## Run
```sh
pip install numpy scipy matplotlib
python model.py            # baseline simulation
python miracle_hunt.py     # search for durable suppression strategies
```

> This is a mathematical model **for study, not medical advice.**

## License
MIT
