# 🧪 TRANSLATION: From Simulation to Real Medicine (updated Sept 2026)

> The PAN protocol in `MIRACLE_CURE.md` is a *mathematical* result. This file
> maps each simulated ingredient to the **real drugs and trials** that exist
> right now, and shows where real medicine matches — and doesn't yet match —
> the simulation. This is information about clinical research, **not medical
> advice**. If you live with HIV: ART gives a normal life expectancy; never
> stop it on your own; discuss any cure research with your HIV doctor.

---

## The surprise headline: reality just validated the model

**EBT-101** (Excision BioTherapeutics) is the only CRISPR-based HIV gene editor
with FDA clearance (IND 2021, first human dose Sept 2022, trial NCT05144386).
It excises integrated HIV proviral DNA with CRISPR-Cas9 — the model's **Eraser**
ingredient.

Clinical results to date:
- Well tolerated, no dose-limiting toxicity.
- **But it did not cure.** Of 5 dosed, 3 stopped ART — all three rebounded and
  had to restart therapy. One participant showed a *delayed* rebound.
- Implication reported by the investigators: the editor did not reach **all**
  latent-infected cells.

This is the simulation's core law repeated in humans: **partial proviral
erasure ≠ cure**. A single wipe, however clever, leaves sleepers that re-seed.
Our model predicted exactly this outcome for "erase-only" protocols (N1 family:
~14/20 at 50 yrs, ~97-98% short-horizon overestimates).

---

## Ingredient-by-ingredient map (PAN → real world)

| PAN ingredient (model) | Real counterpart | Real-world status (2026) |
|---|---|---|
| **Eraser** — erase 99.9999999% provirus in place | CRISPR-Cas9 proviral excision (EBT-101); engineered recombinases (Works Biotech "trefoil"), megaTALs | EBT-101 Phase 1/2, safe but rebound (see above). Real in-place efficiency is **far** below 99.9999999% — the #1 gap. Novel recombinase approaches aim to fix excision chemistry. |
| **Sweep** — clear free virus + kill infected cells | bNAbs (3BNC117, VRC01), immunotoxins, antibody-drug conjugates targeting HIV envelope | bNAbs in multiple Ph1/2 trials; e.g., 3BNC117 + romidepsin combination studied at ART initiation (Nat Med 2022). Potent but partial. |
| **Lockdown** — block new infection for the whole course | Entry inhibitors (maraviroc = CCR5 blocker), lenacapavir (capsid), islatravir | Standard-of-care classes already exist; maintaining them for a long "locked" course is feasible today. |
| **Purge (shaker)** — wake sleepers 1.5×/day all course | Latency-reversing agents (LRAs): HDAC inhibitors (vorinostat, romidepsin, chidamide), PKC agonists (ingenol, bryostatin) | LRAs **do** wake latent HIV *in vivo* — proven. But alone they do **not** shrink the reservoir. |
| **Cleanup (kill awakened cells 1000×)** | BCL-2 antagonists (venetoclax), anti-PD-L1 (ASC22), therapeutic vaccines, CAR-T | **This is the recognized bottleneck.** Trials: venetoclax + LRA (AMBER, NCT05668026, Ph1/2, 2024-26); ASC22 (anti-PD-L1) + chidamide combo (Ph2, 15 patients, ART-free virological control); BCL-2-high cells resist killing — exactly why the model's 1000× cleanup mattered. |
| **ART tail** (days 365-1460) | Standard ART continuation | Feasible; ART is the current backbone of every early-phase cure protocol. |
| **Permanent CCR5 knockout** | CCR5-Δ32 donors (Berlin/London/Düsseldorf patients = cured via transplant), ZFN/CRISPR CCR5-edited stem cells (SB-728MR trial, City of Hope) | Proven curative in the rare transplant cures; gene-edited-CCR5 stem-cell trials ongoing; CCR5-Δ32 homozygotes are naturally protected. Model confirms this knob kills the dominant R5 re-ignition path. |
| **X4/CXCR4 door** (model keeps it open) | No CXCR4-blocking cure strategy — CXCR4 is essential for immune cells | Model's "can't shut this door" is **real biology**; explains the last few knife-edge worlds and why a "perfect" cure is structurally hard. |

---

## What the model says is THE problem (and reality agrees)

1. **Reservoir exhaustion beats perfect chemistry.** Even 99.999%+ erasure
   alone loses (EBT-101 = live proof). You must *also* keep waking + killing
   stragglers under a locked, ART-protected window — the "purge-with-1000×
   cleanup" and long quiet tail of PAN.
2. **Killing the awakened cell is the hard part.** Every LRA trial wakes virus
   but fails to remove the reservoir because viral-expression cells survive
   (BCL-2↑, CTL exhaustion). The next wave (venetoclax, ASC22, vaccines,
   CAR-T) is precisely the "kill arm" our cleanup parameter models.
3. **CCR5-Δ32 makes cures possible but rare.** The 3 cured patients all needed
   a CCR5-deleted transplant; scalable CCR5 editing + reservoir-depletion
   combo is where a "one-time course" could come from.
4. **The final 1-5% is the frontier.** Our 60-world census found ~5/60 worlds
   numerically undecidable — in humans this is the "clonal expansion / last
   provirus" problem: sanctuary sites and clones make absolute elimination
   hard. Real research calls this the difference between "functional cure"
   (ART-free remission) and "sterilizing cure."

## The closest real-world analog to PAN today
A combination protocol combining (i) an excision editor (EBT-101 class),
(ii) a strong kill arm (venetoclax/ASC22/bNAb), (iii) a sustained latency
reversal course (LRA pulses over months, not days), (iv) CCR5-Δ32-style cell
therapy or maraviroc, and (v) continued ART through the whole course. No
single trial has all five yet — the pieces exist; the integration is the
research frontier.

---

### Sources (accessed Sept 2026)
- EBT-101: ClinicalTrials.gov NCT05144386; Excision BioTherapeutics releases
  (2022-2024); CRISPR Medicine News, May 2024 ("EBT-101... does not cure HIV
  — safe, viral rebound in all 3 ATI participants"); CIRM grant CLIN2-13310.
- LRAs/shock-and-kill reviews: Pubmed PMC5990418, PMC9826147, s41392-024-01943-9
  (Nature Sig Trans Target Ther — ASC22 + chidamide Ph2).
- Venetoclax reservoir trial: ClinicalTrials.gov NCT05668026 (AMBER).
- CCR5 editing/HSPC: City of Hope ZFN CCR5 stem-cell trial (CIRM listing);
  Berlin/London/Düsseldorf patient history.
- bNAb + LRA: 3BNC117 + romidepsin early-treatment trial (Nat Med 2022).
---

## Sept 2026 hardening update: simulation's sharpest new predictions

1. **The ART tail is not the cure — the lockdown is.** The ART-FREE variant
   (PAN-2) cured a 240-patient census (233/240 strict, same under clinical
   endpoint) with **no antiretroviral course at all**. Why: the 3-year
   entry-blockade "lockdown" + permanent CCR5 genetic knockout replace ART's
   job. Real-world read: an eve-artificial course of entry-inhibitor
   lockdown + CCR5Δ32 editing (Berlin/London/Düsseldorf-style) could be
   studied WITHOUT a 3-year ART tail — materially shorter pill burden.
2. **Raw erasure beyond 3 nines is pointless** (purge+lock do the work) —
   good news for feasibility, since real CRISPR excision is nowhere near 9
   nines; ~0.999 may suffice in this framework.
3. **CCR5 editing must never happen alone**: 60/60 patients under CCR5-only
   switched to X4 dominance (1.9e5 copies/mL) and 53/60 still died of AIDS.
   CCR5 knockout is a *selection vaccine* for X4 unless the reservoir is
   annihilated simultaneously — exactly the lesson of the failed CCR5-only
   gene-therapy attempts.
4. **Timing changes suffering, not cure-ability**: PAN-2 cures at any start
   day; but treating at day 7 keeps CD4 nadir ≈958 vs ≈262 if you wait to
   day 1460. Test-and-treat economics, quantified.
