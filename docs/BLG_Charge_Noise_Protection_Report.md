> Converted from the author's Word report `BLG_Charge_Noise_Protection_Report.docx` (English). File names in the index refer to the author's working directory; only the scripts present in `scripts/` are published in this repository.

**Charge-Noise Protection in Bilayer-Graphene**

**Quantum-Dot Qubits**

*What Works, What Doesn’t, and Where the Boundary Lies*

A complete, chronological record of the investigation

4-band k·p solver — full Floquet — symmetry analysis — two-electron CI

**Computational research note**

*All numerical results reproduced on an independent machine*

## Contents

## 0.  How to read this document

This note records an entire investigation in the order it actually happened: from an initial manuscript review, through the construction and validation of a numerical solver, the refutation of the manuscript’s central claim, a long sequence of attempts to find a genuinely new protection mechanism, and finally a clear statement of where the boundary of the present tool lies. Nothing has been omitted or smoothed over.

The guiding discipline throughout was simple: **every quantitative claim must pass an explicit anchor** — a symmetry identity, a known analytic limit, a convergence check, or an independent reproduction — or it is reported as a negative result. This is exactly the discipline that separates a durable finding from a retracted one, and it is why several exciting-looking leads in this document are ultimately reported as artifacts.

**Reader’s shortcut.** The single most important outcome is in Section 8 (the structural principle) and Section 9 (the boundary). Sections 3–7 are the evidence that earns them.

## 1.  Where we started: the manuscript

The investigation began with a draft manuscript (“BLG_QD_AC_Protection”) and the question “is it good enough to publish?” The manuscript proposed that applying an AC (periodic) drive to a bilayer-graphene (BLG) valley qubit protects it against charge noise, extending the inhomogeneous dephasing time T₂* by roughly 9% at a drive amplitude of about 50 meV.

The core references were checked first and found genuine and correctly used: Benito et al. (PRB 100, 125430, 2019), Hajati & Burkard (dynamic sweet spots), and the Hecker-group charge-noise measurements in BLG. So the manuscript was not built on fabricated citations — the most common fatal flaw was absent.

However, two problems surfaced immediately from the manuscript’s own numbers. **First**, the mechanism was described as resting on a *concave* splitting curve (∂²Δ₀₁/∂Δ_ext² < 0) invoking Jensen’s inequality — but the manuscript’s own Table 1 values give a slope that becomes less negative with increasing field, i.e. a *convex* curve at the operating point. **Second**, the effect was small (9%) yet required an enormous drive (A_AC = 50 meV, ~2.5× the dot depth), and the supporting appendices were missing from the file.

**Initial assessment:** competent and methodologically careful, but not publishable as written — the mechanism appeared mislabeled, and the central numbers had no demonstrated computational source in the provided files. Verifying them required building the calculation from scratch.

## 2.  Building and validating the solver

Because the manuscript’s numbers could not be traced to any working code in the uploaded files (the provided notebooks contained only failed atomistic attempts and an unrelated 1-D toy chain), the entire calculation was rebuilt from first principles in pure NumPy/SciPy.

### 2.1  Model

A continuum **4-band McCann–Koshino k·p Hamiltonian** in the basis (A1, B1, A2, B2), expanded in a plane-wave basis (box L_box, momentum cutoff q_max), with a Gaussian scalar confinement V(r) = −U₀·exp(−r²/2L_dot²). Standard parameters: a = 0.246 nm, γ₀ = 3.16 eV, γ₁ = 0.38 eV, γ₃ = 0.38 eV.

### 2.2  Five structural anchors — all passed

| # | Anchor | Result |
|---|---|---|
| 1 | Hermiticity | max\|H − H†\| = 0 (machine precision) |
| 2 | Bulk conduction min vs analytic (γ₃=0) | Δ=0.1 eV: numeric 45.5 meV vs analytic 48.4 meV |
| 3 | Trigonal-warping C₃ | dominant angular harmonic m = 3; isotropic when γ₃→0 |
| 4 | Δ_ext = 0 | gapless (as required) |
| 5 | K vs K′ Kramers | spectra identical (exact) |

### 2.3  Confined-dot anchor: Table 1 reproduced exactly

With the same plane-wave box the manuscript implied (N_pw = 293, Hilbert dimension 1172), the orbital splitting Δ₀₁ versus external field reproduced the reference table to all printed digits:

| Δ_ext (eV) | 0.050 | 0.100 | 0.150 | 0.167 | 0.200 | 0.250 | 0.300 |
|---|---|---|---|---|---|---|---|
| Δ₀₁ (meV) | 5.130 | 3.471 | 1.524 | 0.989 | 0.347 | 0.056 | 0.009 |

**Conclusion of Section 2.** The manuscript’s quantitative bulk physics is real and reproducible — not fabricated. An earlier suspicion that the numbers had “no source” was therefore withdrawn. What remained to test were the mechanism and the AC result.

*Files: *blg_kp.py*, *blg_dot.py*.*

## 3.  The mechanism error: convex, not concave

With a working solver, the curvature at the operating point was measured directly:

**Δ₀₁(0.150, 0.167, 0.184) = 1.524, 0.989, 0.593 meV   →   ∂²Δ₀₁/∂Δ_ext² ≈ +0.48  (CONVEX)**

The operating point is convex, not concave. The manuscript’s “concave → Jensen → protection” narrative had the sign of the curvature backwards. The rise of the period-averaged ⟨Δ₀₁⟩ with drive amplitude (its Table 2) is itself the signature of convexity. The numbers were self-consistent; the explanation attached to them was not.

## 4.  Refuting the AC claim with full Floquet

The manuscript’s positive result rests on a period-averaging approximation: it averages the instantaneous qubit splitting over one drive cycle, implicitly assuming the system stays a clean two-level qubit. To test that assumption we evolved the two qubit states through one full drive period in the complete 1172-dimensional Hilbert space (Strang–Trotter in the eigenbasis, no projection) and measured the leakage out of the qubit subspace.

### 4.1  Anchors

- A = 0 → leakage ≈ 0 and Floquet splitting = static Δ₀₁ exactly (7×10⁻¹⁴ meV).

- Time-step convergence: leakage at A = 30 meV is 0.6636 / 0.6633 / 0.6632 for n_step = 200 / 300 / 500.

### 4.2  Result (converged; reproduced on a second machine)

| A_AC (meV) | Leakage | Status |
|---|---|---|
| 5 | 2.9% | 2-level OK |
| 10 | 11% | marginal |
| 15 | 23% | 2-level BROKEN |
| 20 | 38% | 2-level BROKEN |
| 30 | 66% | 2-level BROKEN |
| 50 | 95% | 2-level BROKEN |

At the amplitudes where protection was claimed (A ≳ 20 meV), the drive ejects 38–95% of the qubit population out of the computational subspace every period. The period-averaging method is invalid precisely there. Where the qubit survives (A ≲ 5 meV) the effect is negligible (≲0.4%). **The central positive claim of the manuscript does not survive a full Floquet treatment.**

*Files: *blg_ac_fast.py*, *blg_floquet_full.py*.*

## 5.  The search for a real protection mechanism

With the AC claim refuted, the goal shifted to the underlying question: can a BLG dot qubit be charge-noise protected at all, and if so, how? We pursued this through five model layers. Each subsection records one attempt, its anchor, and its honest outcome.

### 5.1  Single-dot static sweet spot — none

Scanning the orbital splitting over (Δ_ext, U₀), the charge sensitivity |∂Δ₀₁/∂Δ_ext| never reaches an interior zero at a useful splitting; it only shrinks as Δ₀₁ → 0 (the qubit vanishing). No genuine sweet spot. A single dot has no competing-term structure to create an extremum.

*Files: *blg_sweetspot.py, blg_sweetspot2.py*.*

### 5.2  Double-well detuning sweet spot — real (symmetry)

In a symmetric double well, ∂Δ₀₁/∂ε = 0 at ε = 0 to machine precision (≈10⁻¹⁵), independent of box size. This is the BLG transcription of the known Benito–Burkard flopping-mode sweet spot, and it **survives** because it is forced by mirror symmetry. Crucially, it protects only against *detuning* noise — the global displacement field D_ext remains unprotected (slope ≈0.031, marginally worse than the single dot).

*Files: *blg_doubledot.py, blg_compare.py*.*

### 5.3  The d ≈ 60 nm illusion — box artifact

A double-well configuration near d = 60 nm appeared to reduce the D_ext sensitivity ~14×. Box-convergence testing destroyed it: the splitting itself is not converged (Δ₀₁ = 0.99 → 0.31 → 0.14 → 0.16 meV as L_box grows), so neither is the slope. The weakly-bound, widely-separated states leak to the boundary. First clear instance of the pattern that would repeat.

*Files: *blg_dext_scan.py, blg_d60.py, blg_d60_verify.py, blg_converge.py*.*

### 5.4  Valley orbital magnetic moment — method error caught

Turning to the valley qubit, whose splitting derives from the topological orbital magnetic moment μ_orb (valley Zeeman), an initial calculation suggested ~3× lower charge sensitivity than the orbital qubit, and even a sweet spot. Both were wrong: μ_orb had been computed with the unbounded real-space position operator (x·J), which is box-divergent (μ_orb drifted from −2.98 to −7.38 as the box grew). The TRS anchor μ(K) = −μ(K′) held, but the magnitude was meaningless.

The correct gauge-invariant velocity-sum formula gave μ_orb ≈ 0.4 (TRS-exact), an order of magnitude smaller, and still ~20–27% box-drifting. The apparent valley sweet spot’s extremum location **ran away with box size** (0.11 → 0.135 → ≥0.16), never converging.

*Files: *blg_valley.py, blg_moment_correct.py, blg_moment_scan.py*.*

### 5.5  Spin-valley qubit — a real symmetry cancellation

Embedding the spin-valley structure (Kane–Mele λ·ξ·s_z, valley Zeeman via μ_orb, spin Zeeman) in a 4-level model, the fractional charge sensitivity of three qubit types came out as:

| Qubit type | Fractional D_ext sensitivity (/eV) | Protected? |
|---|---|---|
| Orbital | 27.9 | no |
| Valley (Δξ) | 8.5 | partial |
| Spin-valley | 8.2 | partial |
| Spin (Δs, same valley) | 0.00 (machine zero, all fields) | YES (symmetry) |

The spin qubit within a fixed valley is first-order immune to charge noise: the (D_ext-dependent) valley Zeeman shifts both spin states equally and cancels in the splitting. This **survives** because it is a valley-diagonal symmetry. It is, however, an expected consequence (spin does not couple to charge directly), and its residual sensitivity enters only through any field-dependence of the spin–orbit coupling — negligible for intrinsic Kane–Mele in pristine BLG, non-zero for proximity-enhanced SOC (measured to tune 40→80 μeV with field).

*Files: *blg_spinvalley.py*.*

### 5.6  Two-electron sector — a fundamental model limitation

Opening the two-electron layer to seek a singlet-triplet qubit, the envelope-model Coulomb exchange between the two valleys came out as 56 meV — about 160× the measured value (~0.35 meV). The reason is physical and important: the intervalley exchange requires momentum transfer ≈ |K−K′| ≈ 17 nm⁻¹, where the Coulomb kernel 2πe²/ε|q| is small. The envelope model factors out the fast Bloch phase and therefore cannot represent this momentum — it overestimates the exchange catastrophically.

Correcting the kernel to the true |K−K′| momentum gives J_iv ≈ 0.2 meV, matching the measured 0.15–0.5 meV — confirming the missing physics is real, large, and outside the envelope model. The same-valley orbital exchange that IS computable is tiny (~1 μeV) and at the numerical floor.

A two-electron sweet spot (extremum of J_iv vs D_ext) appeared at one box size and vanished at larger boxes — the sixth such artifact.

*Files: *blg_2e.py, blg_2e_v2.py, blg_intervalley.py, blg_iv_converge.py*.*

## 6.  Six candidate sweet spots, all artifacts

Collected in one place, every candidate that relied on tuning a parameter to an extremum failed its convergence check:

| # | Candidate | Apparent effect | Killed by |
|---|---|---|---|
| 1 | AC-driven orbital qubit | ~9% T₂* | full Floquet (38–95% leakage) |
| 2 | Double-well d≈60 nm | ~14× | box size (Δ₀₁ not converged) |
| 3 | Valley moment 3× | ~3× | unbounded operator; box-divergent |
| 4 | Valley charge sweet spot | ~10²× | extremum runs with box |
| 5 | Two-electron S-T sweet spot | sweet spot | vanishes at larger box |
| 6 | (corrected μ_orb extremum) | — | location not box-convergent |

## 7.  Trying new symmetry axes

Having learned that only symmetry survives, the strategy changed: instead of hunting another tuned extremum, look for a symmetry not yet exploited. Three were tried.

### 7.1  Particle-hole symmetry — approximate, and needs p–n geometry

BLG has an approximate particle-hole symmetry, but the γ₁ dimer term breaks it (electron/hole bound-state energy ratio ≈ 0.85, i.e. ~15% breaking). A genuine electron–hole qubit also needs a p–n confinement geometry, not a single scalar well. Approximate symmetry yields only approximate protection.

*File: *blg_ph.py*.*

### 7.2  Intervalley scattering — outside the solver

The valley qubit’s protection is limited by intervalley scattering, a lattice-scale process requiring momentum ≈ |K−K′|. Our basis expands around each valley with |q| < 0.5 nm⁻¹ and cannot represent |K−K′| ≈ 17 nm⁻¹ — the same wall as the two-electron exchange. Improving valley protection is therefore not computable in this tool.

### 7.3  Minivalley C₃ triplet — basis breaks C₃

The solver does reproduce the valley-center → minivalley transition (threefold near-degeneracy at strong gap / large dot, matching Mayer–Knothe). But the square plane-wave grid does not respect C₃ (the bare triplet is already split by ≈0.007 meV — a basis artifact), so C₃-based protection cannot be cleanly studied here. Taken at face value, a central charge fluctuation splits the triplet roughly linearly anyway.

*Files: *blg_minivalley.py, blg_minivalley_protect.py*.*

## 8.  The structural principle

**In electrostatically confined bilayer-graphene quantum dots, within the envelope k·p model, charge-noise protection is** **symmetry-locked****. It comes only from exact symmetries (mirror parity of detuning; the valley-diagonal cancellation of the valley Zeeman for spin qubits; TRS equality of the two valley densities). Tunable “magic operating points” — sweet spots from accidental extrema — are** **not robust****: they are artifacts of finite numerical domains and do not survive convergence.**

This is supported by nine independent attempts: six tuned extrema, all of which dissolved under convergence (Section 6), and three exact-symmetry features, all of which survived (Sections 5.2, 5.5, plus the TRS valley-density equality). The practical corollary for experiment and theory: do not chase tuned charge sweet spots in BLG dots — rely on and engineer symmetry.

## 9.  The boundary, and the concrete next step

“How do we improve the symmetry protection?” is the right question, and its answer defines the boundary of this work. The protection of each symmetry-protected qubit is limited by the term that breaks its symmetry — and all three such terms live outside the present solver:

- Valley qubit → limited by intervalley scattering (momentum |K−K′|, outside small-q basis).

- Particle-hole qubit → symmetry only approximate (γ₁), and needs p–n geometry.

- Minivalley C₃ qubit → square grid breaks C₃; needs a C₃-respecting basis.

All three roads terminate at the same place: **intervalley and C₃ physics that an envelope / small-q / square-grid solver cannot represent.**

### 9.1  Concrete next step

Improving charge-noise protection requires a full-Brillouin-zone atomistic (or full-BZ momentum-space) calculation that represents both valleys at their true separation and respects C₃ — e.g. the dense BZ sampling used in recent work treating BLG dots as Haldane spin-1 simulators (~901×901 k-points, equivalent to ~3×10⁶ atoms). With that tool one can compute, as functions of dot design (size, gap, shape, minivalley vs valley-center regime):

- the intervalley scattering matrix element that limits the valley qubit, and whether the minivalley regime suppresses it;

- the C₃-triplet response to realistic local charge noise;

- the two-electron singlet–triplet structure with correct intervalley exchange.

These are the quantities in which a genuine, symmetry-rooted improvement could live. They are beyond the present envelope solver, but they are well-posed, and the present results define exactly where to begin.

## 10.  File index

| File | Purpose |
|---|---|
| blg_kp.py | 4-band k·p Hamiltonian + five structural anchors |
| blg_dot.py | plane-wave confined dot; reproduces reference Δ₀₁ table |
| blg_ac_fast.py | period-average AC analysis (the flawed method) |
| blg_floquet_full.py | full-Hilbert Floquet + leakage (refutes AC claim) |
| blg_doubledot.py / blg_compare.py | double well; detuning sweet spot (survives) |
| blg_dext_scan / d60 / d60_verify / converge | d≈60 nm artifact and its convergence death |
| blg_sweetspot.py / blg_sweetspot2.py | single-dot sweet-spot scan (none found) |
| blg_valley.py | orbital magnetic moment (TRS anchor) |
| blg_moment_correct.py / _scan / _scan2 | corrected gauge-invariant moment; convergence |
| blg_spinvalley.py | spin/valley/orbital sensitivity; spin cancellation (survives) |
| blg_2e.py / blg_2e_v2.py | two-electron Coulomb matrix elements |
| blg_intervalley.py / blg_iv_converge.py | corrected intervalley exchange; artifact check |
| blg_ph.py | particle-hole symmetry check (approximate) |
| blg_minivalley.py / _protect.py | minivalley regime (reproduced); C₃ basis limit |
| RESEARCH_NOTE.md | machine-readable summary of all of the above |

**Conventions.** *All energies in meV, lengths in nm. Effective dielectric ε_r = 4 where Coulomb enters. The model is spinless / effective single-particle k·p with quasi-static noise unless stated; spin–valley structure (5.5) is added phenomenologically with literature parameters (λ_KM ≈ 30 μeV, g_s = 2, g_v calibrated to ~30). Every numerical result above was reproduced on an independent machine.*
