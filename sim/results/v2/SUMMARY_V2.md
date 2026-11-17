# Glint simulation, round 2: verdict

**How to read this page.** All numbers are from simulation, at 780 nm, using the paper's measured reflectance (see Task 6), with 10 mW/cm² of light. "Read" means the decoder recovers a chip's 8-bit ID and 8-bit value (Task 4). "m" is modulation depth: how much of the light the camera sees blinks. Full tables are in `tables_v2.md`; the data are in the CSVs.

## 1. The whole-plate claim

**Claim tested:** "one ordinary camera reads hundreds of chips across a 96- or 384-well plate at once." **As stated, it is not supported.**

**What does work: one fixed camera on flat chips.** We modelled one 20 MP sensor behind a large object-side telecentric lens (110 × 75 mm field, NA 0.02, 20.6 µm pixels, 20 fps, chips blinking at 5 Hz). It reads every chip on the plate in parallel:
- chips on top of or under organoids: 12 s per read (m ≈ 0.26)
- 100 µm SPOTs buried 300 µm deep: 50 s

96 and 384 wells take the same time, because the plate footprint is the same.

**What breaks it: tilt.** A mirror SPOT tilted by more than about **0.5°** sends its reflection outside an NA 0.02 lens (about 2° at NA 0.05). With random orientation, **0%** of surface chips are readable. If chips settle flat on the glass, half land face-down (readable from below) and half face-up. Within that half, only those flat to ±0.5° work at NA 0.02, which is 13% of all chips. At NA 0.05 the tolerance is 2°, so all the face-down ones work (50% of chips). Also, the telecentric lens is a large special optic, not an "ordinary camera" lens.

**The cheapest robust option is a scanning head (NA 0.25).** It reads a chip in 0.41 s and covers a 96-well plate in about 90 s (384 wells in about 6 min), with 0.5 s per move assumed. It tolerates about 10° of surface-chip tilt.

**Middle ground: 4 or 16 fixed cameras.** These raise the NA to 0.04 or 0.08, which relaxes the tilt limit to about 1.7° or 3.5° and keeps 12 s reads.

## 2. Can a chip at the centre of an organoid be read?

Organoids are modelled as spheres; times are per read.

| organoid diameter | scanning head (NA 0.25), from above or below | whole-plate camera (NA 0.02) |
|---|---|---|
| 0.5 mm | yes, every SPOT size and density (≤ 0.4 s) | loose tissue: yes. Average tissue: only SPOT ≥ 100 µm. Dense tissue: no |
| 1.0 mm | yes, every SPOT size and density (≤ 49 s) | loose tissue only |
| 2.0 mm | yes for a 300 µm SPOT (≤ 6 s) and a 100 µm SPOT in loose or average tissue (≤ 43 s); 30 µm only in loose tissue | no |

**Above vs below.** Reading from below gives about half the modulation, because the dish bottom reflects 4% of the light. It still works.

**Slab vs sphere.** Compared with the sphere, the round-1 slab:
- **under**estimated contrast at the centre by up to **1.7×**, for small, dense organoids (range 0.93–1.67×)
- **over**estimated the absolute signal by up to **3.6×** for large, dense ones

## 3. Chip orientation, and the best physical layout

**Mirror SPOT, no tissue in the path, NA 0.25.** Signal falls to 43% at 5°, 3% at 10°, and zero by 20°. Only **0.9%** of randomly oriented chips stay readable.

**Mirror SPOT embedded ≥ 50 µm deep.** The light reaching it is already diffuse, so tilt costs signal but the scanning head still decodes every orientation. A plate camera decodes almost none (≤ 0.5%).

**Lambertian (diffuse) SPOT, NA 0.25.** Contrast is 5× lower (m 0.053 vs 0.25), but **16%** of random orientations are readable. At NA 0.02 it is unreadable. That would be a new device; it is not demonstrated.

**Best layout: chip resting on the dish, SPOT facing down into the glass, read from below, sensor side facing up into the organoid.** No tissue sits in the optical path.
- m = **0.26**
- read time 0.41 s with the scanning head, 12 s with a plate camera
- if the chip is flat to ±0.5° (NA 0.02) or ±2° (NA 0.05)

What limits it:
- the 4% dish-bottom reflection, the largest background term
- how flat the chip lies
- the SPOT still needs electrolyte between it and the glass, which a spacer must provide (not modelled)
- focus is easy (depth of field is millimetres at NA 0.02)

## 4. Decoder and the measured noise floor

We built a decoder and tested it on synthetic camera frames. Each chip sends a 29-bit code, and the frames include all the requested noise sources. Using a reference region to cancel light-source flicker, the lowest m that decodes reliably (≥ 23/24 chips) with the scanning head was:
- 6×10⁻⁴ to 1.2×10⁻³ in a 2 s read
- 2.5–3×10⁻⁴ in 10 s
- about 1.3×10⁻⁴ in 60 s

These floors are the same at 0.01%, 0.1% and 1% source noise. Without the reference region, 1% noise raises them about 20–40×. The whole-plate camera's floor is 1.2×10⁻² (12 s) or 5×10⁻³ (56 s), because it collects far fewer photons per second.

**The floor is not one number.** It is set by photons per read: "SNR per blink cycle × √(cycles per bit) ≥ 23", or ≥ 120 for chips deeper than 500 µm, where this decoder is less efficient. That rule agrees with all 72 decoder runs.

**Updated maximum depth** (a chip decodes within 60 s, scanning head, 780 nm):

| SPOT | loose tissue | average tissue | dense tissue |
|---|---|---|---|
| 30 µm | ≥ 1000 µm | 500 µm | 500 µm |
| 100 µm | ≥ 1000 µm | ≥ 1000 µm | 500 µm |
| 300 µm | ≥ 1000 µm | ≥ 1000 µm | ≥ 1000 µm |

With a flat-chip whole-plate camera, the limit is 100–300 µm in average or dense tissue (500 µm in loose). The round-1 assumed floor (10⁻³) gave 300–500 µm. The decoded floor is lower only because our noise model contains nothing near the blink frequency except shot noise. **Real biological and electrical noise could raise it.** The bench experiment measures this.

## 5. Light dose

During reads the light is 10 mW/cm².

| how it is read | every 1 min | every 5 min | every 15 min |
|---|---|---|---|
| scanning head, only the well being read is lit (0.46 s) | 6.7 J/cm²/day (0.77% duty) | 1.3 J/cm²/day | 0.45 J/cm²/day |
| fixed plate camera (12 s reads) | 179 J/cm²/day | 36 J/cm²/day | 12 J/cm²/day |
| scanning head, whole plate lit during the 90 s scan | 864 J/cm²/day | 253 J/cm²/day | 84 J/cm²/day |

**What the one source we could read found.** Khan, Tang & Arany (Sci. Rep. 5:10581, 2015) saw 810 nm phototoxicity in keratinocytes at ≥ 27 J/cm² delivered at 90 mW/cm². That happened only in light-absorbing black wells, and they tied it to heating to ≥ 45 °C plus reactive oxygen; clear wells showed no significant toxicity.

**Why we make no safety judgement.** That study is a single, more intense exposure, not chronic daily low-level light. We found no threshold for our regime. We could not access the photobiomodulation dose literature (Huang et al., 2009 review), which reports biological effects of red/NIR light at low doses. So we report the doses and make no safety judgement.

**Energy storage.** None is needed. The chip is lit during its own read, which is when it needs power: about 1 µW on a 100 µm cell, a 50 ms start-up assumed. Storage would matter only for sensing or logging between reads.

## 6. Which ΔR, and from where

We used R_on = 0.463 and R_off = 0.342, so ΔR = 0.121. These come from the paper's sigmoid fit:
- **Paper:** arXiv:2606.23916, Suppl. Section S4, Fig. S11, Eq. S40
- **Device:** 45 nm PEDOT:PSS on Pt, 50 µm SPOT, 1× PBS
- **Measured at:** 635 nm
- **Read at:** the +200 / −800 mV drive

The paper gives no number at 780 nm, and says 635 nm is near the wavelength of largest contrast. So the 780 nm headline is **optimistic**, and this is the most important number to measure (bench plan). This ΔR sits inside round 1's 0.05–0.2 sweep. The lower mean R (0.40 vs 0.50) raises m.

## What we can honestly claim

"A validated light-transport simulation supports the feasibility, under stated assumptions, of reading SPOT chips placed on or under an organoid: the blink changes the light a camera sees by about 26%, and a decoder recovers a chip's ID and reading in under a second with a scanning reader. It also supports reading chips at the centre of organoids up to 1 mm across with that reader, and, for chips that lie flat, reading a whole plate in parallel with one camera."

## What we cannot claim

- That the device works. A simulation supports feasibility under stated assumptions; it does not prove the device works.
- That **one ordinary camera** reads a whole plate. That needs a large telecentric lens and chips flat to about 0.5°, or else several cameras or a scanner.
- That randomly oriented chips work with a plate camera.
- Any ΔR at 780 or 850 nm.
- Real organoid optics.
- That deep chips read in real life. The decoded floor contains no biological or chip-electrical noise.
- Cell safety of the light dose.
- That a chip can sit SPOT-down on glass and still switch (electrolyte access).
- That a Lambertian SPOT is possible.

## Three sentences for a pitch slide

1. In our simulation, a chip sitting on or under a mini-organ makes the light a camera sees blink by about **26%**, enough to read its ID and measurement in **under half a second**.
2. Even at the **centre of a 1-millimetre organoid**, a 100-micron chip is decoded in **under one second** in simulation.
3. A scanning reader covers a **96-well plate in about 90 seconds**; read every 5 minutes, each organoid is lit for **0.15% of the day**.
