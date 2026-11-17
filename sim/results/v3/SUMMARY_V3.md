# Glint simulation, round 3 (final): verdict

This supersedes round 2's claims. All results are from simulation.
- **Reflectance.** "635 nm" uses the reflectance change the paper actually measured (R_on 0.463 / R_off 0.342). "780 nm" extrapolates that value and is optimistic.
- **Settings.** Average tissue (μs′ = 1/mm at 780 nm, so 1.3/mm at 635), 10 mW/cm², and the decoder rule calibrated in round 2. "Read" means the ID and value are decoded.
- **Where the details are.** Tables are in `tables_v3.md`.

**Bottom line**
- At the measured wavelength, everything that worked in round 2 still works for chips lying flat on the dish and for flat chips near the surface.
- **Random orientation is the biggest new problem for chips inside an organoid.** Only large SPOTs (300 µm) stay reliable.
- **What should be measured first** is the real noise at the blink frequency inside an organoid, and how chips actually land.

## 1. Headline at 635 nm (measured) vs 780 nm (extrapolated, optimistic)

| result | 635 nm (measured reflectance) | 780 nm (extrapolated, optimistic) |
|---|---|---|
| chip on dish, SPOT down, 100 µm SPOT: scanning head (NA 0.25) / plate camera (NA 0.02, flat) | 0.41 s / 12 s (m = 0.25) | 0.41 s / 12 s |
| whole plate, scanning head: 96 / 384 wells | 88 s / 350 s | 88 s / 350 s |
| 100 µm SPOT 300 µm deep (flat): scanning head / plate camera | 0.41 s / 860 s | 0.41 s / 50 s |
| flat chip at the organoid centre, readable within 60 s (scanning head, from above / below; of 27 size × SPOT × density cases) | 24 / 22 | 24 / 24 |
| same, plate camera (from above) | 7 of 27 | 8 of 27 |
| max decodable depth, scanning head: 30 / 100 / 300 µm SPOT (average tissue) | 500 / ≥1000 / ≥1000 µm | 500 / ≥1000 / ≥1000 µm |
| max decodable depth, plate camera: 30 / 100 / 300 µm SPOT | 150 / 200 / 250 µm | 200 / 300 / 300 µm |

**Why 635 nm looks slightly worse.** At 635 nm, tissue scatters about 1.3× more than at 780 nm. So the plate camera and the deepest chips do slightly worse at the wavelength that was actually measured. The scanning-head headlines are unchanged.

**If the true swing at 780 or 850 nm is smaller.** We tested a swing of 50% or 25% of the 635 nm value:
- At 50%: no headline changes. Chip-on-dish m falls to 0.12, still 0.41 s.
- At 25%: chip-on-dish m = 0.06 (still 0.41 s). The 100 µm SPOT's maximum depth falls from ≥1000 to 500 µm. Readable centre cases drop from 24 to 22.

## 2. Chip at the centre with random orientation

Scanning head, 635 nm, read from above. A mirror SPOT tilted away from the camera sends its reflection elsewhere.

| organoid | 30 µm SPOT: ≤1 s / ≤60 s | 100 µm SPOT: ≤1 s / ≤60 s (median) | 300 µm SPOT: ≤1 s / ≤10 s / ≤60 s |
|---|---|---|---|
| 0.5 mm | 4% / 21% | 22% / 54% (35 s) | 52% / 69% / 92% |
| 1.0 mm | 0% / 12% | 15% / 61% (29 s; 90th percentile 480 s) | 61% / 97% / 100% |
| 2.0 mm | 0% / 0% | 0% / 2% | 3% / 37% / 70% |

**"A chip at the centre of a 1 mm organoid is decoded in under one second" does not survive random orientation** for a 100 µm SPOT: that is true only 15% of the time.

**Corrected sentence:** *"With random orientation, a 300 µm SPOT at the centre of a 1 mm organoid is decoded within 10 s 97% of the time; a 100 µm SPOT is decoded within a minute only 61% of the time."*

Loose tissue makes orientation matter more, not less, because the light stays directional.

## 3. Flatness and landing, for a chip on the dish with its SPOT facing the glass

**Debris tolerance by reader.** "Largest particle" is the height that keeps the chip readable; the 5 µm row shows the tilt a typical speck causes.

| reader | tilt limit | largest particle, 150 / 200 / 300 µm chip | a 5 µm particle under a 150 µm chip (1.9° tilt) |
|---|---|---|---|
| one plate camera, NA 0.02 | 0.52° | 1.4 / 1.8 / 2.7 µm | **fails** |
| plate camera, NA 0.05 | 2.0° | 5.3 / 7.1 / 10.6 µm | OK, barely |
| 16 cameras, NA 0.08 | 3.1° | 8.1 / 10.8 / 16.2 µm | OK |
| scanning head, NA 0.25 | 10.7° | 28 / 38 / 56 µm | OK |

**Share of settled chips readable.** This assumes tilt spread evenly from 0 to 2° and half of chips landing SPOT-up.

| reader | one SPOT, read from below | (a) SPOT on both faces | (b) read from below and above (chip beside the organoid / under a 1 mm organoid) |
|---|---|---|---|
| NA 0.02 | 13% | 26% | 26% / 13% |
| NA 0.05 or 0.08 | 50% | 100% | 100% / not simulated |
| scanning head | 50% | 100% | 100% / 100% |

**Fix (c), described only.** Make the SPOT face denser than the back (for example a thicker metal layer under the SPOT), or add a thin magnetic layer and put a magnet under the plate. Either biases which face settles down. Whether a 150–300 µm plate sedimenting in medium actually lands face-down this way must be tested.

**Ordinary lab conditions** (a few µm of debris): the scanning head tolerates them, and so do several cameras at NA ≥ 0.08. A single fixed NA 0.02 camera does not.

## 4. What to measure first

Sensitivity ranking: scanning head, 635 nm, 100 µm SPOT. Each input is varied across its plausible range.

| rank | input (range tested) | effect | bench measurement |
|---|---|---|---|
| 1 | added noise at the blink frequency (0 to 10⁻² /√Hz) | centre of 1 mm: 0.41 s → 2.6 h. 2 mm: 44 s → never. Max depth: ≥1000 → 150 µm | lock-in recording of an unpowered chip inside a real organoid |
| 2 | chip orientation (flat → random; dish tilt 0–10°) | centre of 1 mm, 90th percentile: 0.41 → 480 s. Dish: no change up to 10° | image how settled and embedded chips actually lie; wedge test |
| 3 | tissue scattering μs′ (0.3–2.0 /mm) | centre of 2 mm: 0.41 → 320 s. Max depth: ≥1000 → 500 µm | organoid μs′ (OCT attenuation or integrating sphere) |
| 4 | reflectance change ΔR (100% → 25%) | centre of 2 mm: 44 → 710 s. Max depth: ≥1000 → 500 µm | bare-chip reflectance at the reader's wavelength |
| 5 | light level (1–100 mW/cm²) | centre of 2 mm: 44–61 s; nothing else changes (the camera saturates). It still sets the chip's power (about 10 mW/cm²) | power meter at the sample |
| 6 | absorption μa (0.01–0.05 /mm) | no change in any case run | spectrophotometry of an organoid homogenate |

## 5. Pitch assets

- `chart_modulation_vs_depth.png` / `.svg`: 635 nm, average tissue, the decoder floor drawn, and the caption "Simulation, not measurement".
- `decoder_demo.gif` and `decoder_demo_strip.png`: labelled "Synthetic data". The decoder recovered all 3 chips' IDs and values.

## What we can honestly claim (supersedes round 2)

"A validated light-transport simulation supports the feasibility, under stated assumptions, of reading SPOT chips that lie on the dish with the SPOT facing the glass, using the reflectance change measured by the inventors: about 25% blink contrast, decoded in under half a second by a scanning reader that covers a 96-well plate in about 90 s and tolerates a 10° tilt. Inside organoids, it supports reading large (300 µm) SPOTs at any orientation at the centre of a 1 mm organoid, usually within seconds."

## What we cannot claim

- That the device works. A simulation supports feasibility under stated assumptions; it does not prove the device works.
- That **one ordinary camera** reads a whole plate. A single fixed camera needs a large telecentric lens and chips flat to 0.5°, which a few µm of debris breaks.
- That **small (≤100 µm) SPOTs** inside organoids read reliably at random orientation, or that any chip reads at the centre of a 2 mm organoid.
- Any reflectance change at 780 or 850 nm, or real organoid optics.
- That deep chips read in real tissue. The decoded floor contains no biological or electrical noise at the blink frequency, which is the top-ranked unknown.
- That chips land SPOT-down, or still switch while pressed to glass (electrolyte access).
- Cell safety of the light dose.

## Three pitch sentences (635 nm, measured reflectance; orientation as stated)

1. In our simulation, a chip lying SPOT-down on the dish blinks with **25% contrast** and is decoded in **0.4 seconds** by a scanning reader, even when tilted **up to 10°** by debris beneath it.
2. That reader covers a **96-well plate in about 90 seconds**.
3. A **300-micron** SPOT at the centre of a 1-millimetre organoid, at **any orientation**, is decoded within **10 seconds 97% of the time** in simulation.
