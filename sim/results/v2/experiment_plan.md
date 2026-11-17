# Smallest bench experiment that tests the simulation

**Goal.** Test the model's central prediction with a few real SPOT chips: how the blink's modulation depth m and the image blur change as a chip is buried under a tissue-like layer of known scattering. The test uses a microscope-like optic and a whole-plate-like optic. All predictions below come from `sim/bench_prediction.py` and are saved in `bench_prediction.csv`. They model this exact geometry: a chip on glass under a phantom of thickness d.

## Materials

- **Chips.** 2–3 SPOT temperature chips from the inventors, plus, if available, one bare SPOT test die with external drive.
- **Phantom.** Agarose + Intralipid in 1× PBS, so the phantom itself is the electrolyte the SPOT needs.
  - **Recipe** for 10 mL at μs′ ≈ 1.0 /mm (780 nm):
    1. Dissolve 0.10 g agarose (1% w/v) in 9.5 mL 1× PBS by heating.
    2. Cool to about 50 °C.
    3. Add 0.49 mL Intralipid-20% (4.9% v/v), mix gently, and cast.
  - **Basis for the recipe.** The van Staveren et al. (Appl. Opt. 1991) fits, as summarised by Jacques at omlc.org/spectra/intralipid, give Intralipid-10% μs ≈ 2.54×10⁹·λ^−2.4 cm⁻¹ and g ≈ 1.1 − 0.58×10⁻³·λ (λ in nm). That is μs′ ≈ 10.3 /mm at 780 nm for the 10% stock. So about 9.7% v/v of Intralipid-10%, or 4.9% v/v of Intralipid-20%, gives about 1.0 /mm. This assumes scattering scales linearly with lipid concentration.
  - **Measure the batch.** Jacques warns that batches vary. Measure collimated transmission of a 1:50 dilution in a 1 cm cuvette at 780 nm to get μs, then μs′ = μs(1 − g) with g ≈ 0.65. Optional extra batches: 0.5 and 2.0 /mm (half and double the Intralipid).
- **Thicknesses.** 0 (PBS only), 100, 200, 300, 500 and 1000 µm. Cast sheets between glass slides with precision shims or stacked #1.5 coverslips. Measure every sheet with a micrometer.
- **Assembly.** Put the chip on a glass slide, add a PBS film, and lay the phantom sheet on top. Keep it wet.

## Optics

1. **Microscope-like.** 780 nm LED epi-illumination (for example the paper's CHROLIS-type 780 nm channel) at about 10 mW/cm². Use a 10× NA 0.25 objective and an sCMOS camera in ROI mode at ≥ 600 fps. The chip blinks at about 150 Hz.
2. **Whole-plate-like.** An object-side telecentric lens with NA about 0.02 and coaxial 780 nm illumination, on a 20 MP CMOS at about 20 fps. The chip blinks at about 5 Hz, or slower if its firmware is fixed.

## What to measure

1. **Bare chip in PBS at 780 nm.** Measure m and R_on/R_off. This is the single biggest unknown: the paper gives ΔR only at 635 nm.
2. **m under each phantom thickness.** Lock-in at the blink frequency, using a reference region for common-mode correction (as in `decoder.py`). Compute m(d)/m(0).
3. **Blur.** Measure the image's 50%-energy radius (r50) for each thickness.
4. **Noise floor.** Record with the chip unpowered and compute the apparent m at the blink frequency. This replaces the simulation's noise model with a measured one.
5. **Tilt, low-NA optic, no phantom.** Put the chip on wedges of 0.25°, 0.5°, 1°, 2° and 5°.
6. **Decode.** If the chip's code is known, run `decoder.py` on the recorded frames and count correct readouts.

## Model predictions

For a 100 µm SPOT. The range in brackets is μs′ ± 20%. A 50 µm SPOT gives similar values out to 300 µm and falls faster beyond it.

| optic | 100 µm | 200 µm | 300 µm | 500 µm | 1000 µm |
|---|---|---|---|---|---|
| NA 0.25, m(d)/m(0) | 0.91 (0.88–0.92) | 0.74 (0.66–0.81) | 0.49 (0.38–0.62) | 0.12 (0.06–0.24) | 0.0027 (0.0014–0.006) |
| NA 0.25, image r50 (sharp = 40 µm) | 40 µm | 41 µm | 47 µm | 94 µm | 450 µm |
| NA 0.02, m(d)/m(0) | 0.71 (0.62–0.80) | 0.26 (0.14–0.42) | 0.051 (0.019–0.14) | 0.0022 (0.0008–0.009) | ~2×10⁻⁵ |

**Tilt, no phantom.** For a mirror SPOT:
- At NA 0.25, signal should fall to about 43% at 5° and about 3% at 10°.
- At NA 0.02, signal should vanish beyond about 0.5–0.9°.

## Confirm or contradict

**Confirms the model** if all of these hold:
- the measured m(d)/m(0) falls inside the bracketed ranges at 200–500 µm, widened by the measured μs′ uncertainty
- r50 stays near the sharp size up to about 300 µm and grows several-fold by 500 µm
- the tilt cut-offs appear at the predicted angles

**Contradicts the model** if any of these hold:
- m at 300 µm (NA 0.25) is below about a quarter of the prediction, or the blur at 500 µm is more than twice the prediction. Either means the model is too optimistic and every embedded-depth claim must be cut.
- The unpowered-chip noise floor at the blink frequency is far above shot noise, for example apparent m above 10⁻³ in 10 s. Then the decoded floors in SUMMARY_V2 are too optimistic, and deep chips must be judged against the measured floor instead.
- Results are much *better* than predicted. That would also mean the model is wrong, more likely in the phantom's μs′ or in g.
