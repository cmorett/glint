# Glint light-transport simulation: verdict

**Question.** Can a camera read a blinking SPOT chip that sits on, under or inside an organoid, and how deep can the chip be?

**Short answer.**
- **On top (A) and underneath (B):** both are easy. In every swept setting, modulation depth is m ≈ 0.07–0.19, and SNR 5 takes one blink cycle (6.7 ms).
- **Embedded (C), photon budget:** counting photons is not the problem at any simulated depth. With shot noise only, every baseline configuration down to 1000 µm reaches SNR 5 in under 60 s. Across all 19,683 configurations, 27 are marginal and none is unreadable.
- **Embedded (C), what actually limits it:** contrast. The blink becomes a tiny ripple on a large background of scattered light. We use m ≥ 10⁻³ as an assumed practical floor (LED drift, organoid motion, camera fixed-pattern noise). It is a judgement call, not a measured number. Against that floor, the limit is a few hundred µm for small SPOTs in average or dense tissue, and about 1 mm for large SPOTs.

## Maximum depth (µm), embedded chip

These are the depths where m stays ≥ 10⁻³. Baseline settings: NA 0.25, 10 mW/cm², ΔR = 0.1, μa = 0.02/mm. "Readable" in the shot-noise sense holds to ≥1000 µm in every cell, and 1000 µm is our deepest simulated depth.

| SPOT | tissue μs′ (at 780 nm) | 635 nm | 780 nm | 850 nm | 780 nm + λ-shift (filter 10³, eff. 10%) | 780 nm + λ-shift (10⁴, 100%) |
|---|---|---|---|---|---|---|
| 30 µm | 0.3 loose | >1000 | >1000 | >1000 | >1000 | >1000 |
| 30 µm | 1.0 baseline | ~490 | ~640 | ~700 | >1000 | >1000 |
| 30 µm | 2.0 dense | ~300 | ~375 | ~400 | ~970 | >1000 |
| 100 µm | 1.0 | ~780 | ~910 | ~950 | >1000 | >1000 |
| 100 µm | 2.0 | ~490 | ~600 | ~640 | >1000 | >1000 |
| 300 µm | 1.0 | >1000 | >1000 | >1000 | >1000 | >1000 |
| 300 µm | 2.0 | ~870 | ~960 | ~980 | >1000 | >1000 |

Loose tissue (0.3) is >1000 for every SPOT size. A stricter floor (10⁻²) would roughly halve these depths. For example, a 100 µm SPOT in baseline tissue at 780 nm has m = 0.079 at 300 µm and 0.018 at 500 µm.

## Matt's points

1. **Is ~100 µm right for a small SPOT? Yes, if "resolve" means seeing the SPOT as a sharp object.** We measured where the lateral spread of reflected light at the tissue surface (2 × r50) grows larger than the SPOT itself:

   | SPOT | baseline tissue | dense tissue |
   |---|---|---|
   | 30 µm | 70–90 µm | 50–60 µm |
   | 100 µm | 130–160 µm | 90–110 µm |
   | 300 µm | 240–290 µm | 170–200 µm |

   Seeing a blurred blink is a weaker requirement than resolving the SPOT. A focused camera still collects a sharp core of weakly scattered light sitting on a halo. That is why detection reaches deeper than resolution does.
2. **How much does a 300 µm SPOT help?** About 3× deeper than a 30 µm SPOT before blur exceeds its size. At 1 mm in baseline tissue it gives about 8× more contrast than a 100 µm SPOT (m 0.0039 vs 0.00051) and about 80× more than a 30 µm SPOT. The tradeoff is that it is larger than most of the chip.
3. **Do the on-top and underneath placements work comfortably?** Yes, with a wide margin. The underneath placement has about 60% more background, mostly from the dish-bottom reflection (4%), but it still gives m ≈ 0.17.
4. **What would wavelength shifting buy? (Hypothetical.)** It multiplies contrast by roughly (conversion efficiency × filter rejection) for deep chips. With 10% efficiency and 10³ rejection, contrast for a 100 µm SPOT at 1 mm rises ~70×, from 5×10⁻⁴ to 0.036. At 780 nm and NA 0.25, that keeps every simulated depth above the 10⁻³ floor except a 30 µm SPOT in dense tissue, which drops below it at about 970 µm. It buys nothing when efficiency × rejection ≤ 1 (for example 1% with 10²). It is only worth building for embedded chips.

## Power check (our rough estimate: 100 nW electrical, 10% PV, so 1 µW optical on a 100×100 µm cell)

Tissue here absorbs little and mostly redistributes light. The downward irradiance at the chip is 0.84–1.04× the surface value down to 1 mm. The cell therefore needs about **10 mW/cm² at any depth up to 300–500 µm, and 11–15 mW/cm² at 1 mm**. 1 mW/cm² is 10× too little. Whether cells tolerate 10–15 mW/cm² of continuous NIR light is not tested here.

## Sanity check vs the paper

We modelled 250 µm of chicken skin (μs′ = 2.45/mm, Jacques skin values) over a 50 µm SPOT at 780 nm:
- The image is blurred: the 50% radius is 66–78 µm, compared with 20 µm for a sharp SPOT.
- The blink is still detectable: m = 0.003–0.02, and the shot-noise SNR is far above 5 per cycle.

This matches the paper's "blurred but visible" observation. It is a qualitative match only, because the paper reports no contrast numbers.

## What we can honestly claim

"A validated Monte Carlo light-transport simulation supports the feasibility, under stated assumptions, of optically reading SPOT chips placed on top of or underneath an organoid: the blink modulates the camera signal by about 17%. It also supports reading chips embedded a few hundred micrometres deep, with about 1% modulation at 300 µm in dense tissue for a 100 µm SPOT, using about 10 mW/cm² of light, which also delivers about 1 µW to a 100 µm solar cell at those depths."

**We cannot claim:**
- that the device works. A simulation supports feasibility under stated assumptions; it does not prove the device works.
- real organoid optics. These are generic soft-tissue values that we did not measure.
- the paper's actual ΔR. We swept 0.05–0.2.
- reliable reading below m ≈ 10⁻³, or embedded reading at ≥500 µm in dense tissue.
- anything about wavelength shifting. No such SPOT exists.
- cell safety at 10–100 mW/cm².
- one camera reading a whole plate at NA ≥ 0.25. Whole-plate optics have NA ≪ 0.1.
- that the chip doesn't harm or perturb the organoid.
- anything about the on-chip sensors.

## Assumptions (source) and simplifications

**Assumptions:**
- n tissue 1.37, medium 1.33, glass dish 1.50; g = 0.9 (prompt).
- μs′ = 0.3 / 1.0 / 2.0 mm⁻¹ at 780 nm (prompt). Scaled to other wavelengths as λ^−1.286 (Jacques 2013, "other soft tissues" b).
- μa = 0.01–0.05/mm, wavelength-independent (prompt; assumption).
- Slab thickness 2 mm (assumption).
- R = 0.5 ± ΔR/2 (prompt; the paper gives no absolute value we could extract).
- Collimated normal illumination at 1 / 10 / 100 mW/cm² (assumption).
- QE 0.6 and 150 Hz blinking (prompt and paper).
- Camera (assumption): 5 µm object pixel, 30k e⁻ full well, 1.6 e⁻ read noise, 300 fps. Brightness is scaled down until no pixel saturates.
- ROI radius chosen to maximise SNR.
- Glare from the outermost flat interface: 2% (water/air) or 4% (air/glass).
- Contrast floors 10⁻³ / 10⁻⁴ (assumption).
- Power need (our own estimate).

**Simplifications:**
- The organoid is a flat, laterally infinite slab, not a sphere.
- The chip neither shadows nor perturbs the light.
- Ideal aberration-free optics, focused on the chip plane.
- Only shot noise, plus camera full well and read noise. There is no source, motion or biological noise.
- In the wavelength-shift case, the shifted light sees the same tissue optics.

## Validation and statistics

All checks pass (`results/validation.txt`):
- MCML benchmark: Rd 0.0977 / 0.0972 vs 0.0973; Tt 0.6607 / 0.6611 vs 0.6610. That is ≤0.4% error, in both absorption schemes.
- The zero-scattering limit matches Beer–Lambert plus Fresnel to <0.001%.
- Energy is conserved to <10⁻⁵.

The run used 3 seeds × 3×10⁵ photons per stage per setting and took 142 s on 15 cores. CSV errors are standard errors across seeds. The median relative error on m is ≤0.3% down to 300 µm, 1% at 500 µm and 2% at 1 mm. The worst case is 11%, for the smallest signals.
