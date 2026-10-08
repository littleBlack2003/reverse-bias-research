# PM6:Y6 100–300 K closure review

## Verdict

The Gaussian Fermi–Dirac thermodynamic closure and thermodynamically signed finite-volume transport can be extended numerically to 100 K for σe=60 meV and σh=74 meV. This does **not** extend the experimental calibration of mobility, establish contacts, or make the low-temperature device prediction quantitative. Keep three statements separate: exact FD integral, calibrated effective mobility over its measured temperature window, and a conditional full-device stress test.

## Minimal model and its scope

1. Use Gaussian DOS, N0=2.4×10²⁰ cm⁻³, independent 60/74 meV widths, exact Fermi–Dirac occupation. The common N0 is the paper's assumed Y6 molecular density for its two-DOS treatment, not an independent measurement of each transport-network lattice spacing.
2. Generalized Einstein relation D/μ=kBT c/(dc/dη). The secant-Einstein SG face closure and the earlier inverse-activity face closure have the same continuum law and exact equilibrium, but differ at finite mesh. Their finite-mesh values need not match. Both must be mesh-refined.
3. Use effective SCLC zero-field mobility μ300=8.4×10⁻⁴/1.3×10⁻⁴ cm² V⁻¹ s⁻¹ and μ(T)/μ300=exp[−(4/9)(σ/kB)²(T⁻²−300⁻²)]. The missing minus sign in source SI Eq. S1c is a source-level typo, not a MinerU error. This is not automatically the microscopic zero-density limit. Below the measured SCLC window (~223 K), label it extrapolation.
4. Do not append uncalibrated EGDM density and field factors. At 100 K, σe/kBT=6.963 and σh/kBT=8.587, beyond the stated 2–6 temperature window of that empirical transport interpolation. Using N0⁻¹/³≈1.609 nm as hopping distance would add a lattice identification assumption; it is not another experimentally fitted parameter.
5. Use R=β(T)np[1−exp(−A)], A=Δμ/kBT, for reversible dark recombination. This guarantees R=0 at equal quasi-Fermi levels and RA≥0. Prescribed G is a nonequilibrium external light source. The electrical/chemical work identity is testable, but is not a complete photon-energy or total entropy ledger. A fully reversible optical channel needs absorption/emission rates and photon occupation or separately justified photon affinity; these are not supplied by G alone.
6. Contacts, dielectric constant, collection efficiency, and absolute temperature-dependent generation remain separate inputs. A zero DOS-center barrier means occupation 1/2 at the majority contact: an idealized high-injection limit, not a measured work function. The default εr=3.5 is a diagnostic scenario, not a paper-derived calibration.

## Recombination calibration and holdouts

The paper independently reports k2≈8×10⁻¹² cm³ s⁻¹ at room temperature. With constant G and PIA carrier density nPIA(T), β(T)/β300=[nPIA(300)/nPIA(T)]² is a transparent inference. The absolute scale G300=β300 nPIA(300)² remains conditional on matching sample/illumination and on the n=p interpretation.

PIA n(T) is then a calibration input, not a validation result. Reproducing it by construction does not demonstrate a correct recombination model. Likewise, using a bandgap already fitted to Figure 5a makes Figure 5a a reproduction/consistency comparison. Table S3 distinguishes fitted gap from measured disorder; the OCR omitted these footnotes and they must be retained in interpretation. A single fixed exact-FD gap must not be changed at every temperature to absorb discrepancies. The separate S16 lower-intensity voltage curves and separately conditioned BACE data are useful holdouts, with the illumination/device matching caveat.

A fixed Langevin reduction factor is especially unsafe. Corrected zero-field GDM gives [μe(100)+μh(100)]/[μe(300)+μh(300)]=4.166×10⁻⁹, so a fixed-G, fixed-reduction Langevin model would force n100/n300≈15,492. The PIA inference in the present inputs gives only ≈5.90. For the εr=3.5 diagnostic, βPIA/βLangevin is 0.01595 at 300 K, 0.0514 at 200 K, 1.39 at 150 K and 1.10×10⁵ at 100 K. This is a warning about jointly interpreting the low-T effective mobility extrapolation and recombination, not proof that the PIA measurement is wrong or that a universal Langevin upper bound applies to the blend.

## Numerical verification

The FD implementation uses 801-point Gauss–Legendre quadrature of the standard normal over ±16, a dense negative-η half-table with consistent Hermite derivative, particle–hole symmetry, stable logarithmic tails, and monotone inverse seeds with Newton correction. Direct adaptive QUADPACK integration over ±20 is the independent numerical reference.

- 117 EOS checks: 100,120,150,200,223,250,300 K for both widths, zero-disorder limit and maximum validated reduced width 8.7
- Maximum relative density error: 1.16×10⁻¹²
- Maximum relative dc/dη error: 4.41×10⁻¹⁰
- Maximum inverse-density error against independent integration: 1.03×10⁻¹²
- Positive susceptibility, c(-η)=1−c(η), correct dilute and single-level limits, and derivative consistency checked
- 80 independent face-closure checks: exact equilibrium, nonnegative flux–affinity product, and second-order convergence for both SG and inverse-activity closures
- 54 device checks and six solved dark equilibria: current/reaction zero at equilibrium, discrete work telescope, reaction/transport signs, continuum flux, colored Jacobian, and Poisson charge balance

Test code and detailed results are in this directory. `fd_eos_results.json`, `flux_comparison_results.json` and `device_results.json` include numerical tolerances; device results include source hashes. Saved-state and mesh checks must be refreshed against final source/output before packaging.

## Current resolution, especially at 100 K

An absolute current residual such as 10⁻¹² A cm⁻², or a residual relative only to the gross generation/recombination current, cannot certify a much smaller extracted current. Early absolute-quasi-Fermi trial states had wrong-signed dark currents whose spread exceeded their mean. These were numerical zeros, not physical predictions. Contact-referenced quasi-Fermi variables remove a major cancellation source.

The ledger now reports current spread relative to mean current, a closure-based uncertainty indicator (spread + integrated-continuity error + a floating-point cancellation allowance), and whether the sign is resolved. A 0.1% relative-current target with an explicit machine-roundoff allowance is included in the numerical gate. The indicator does not include mesh, parameter, or material uncertainty. At open circuit, convert the absolute current uncertainty to a voltage uncertainty through a locally measured dJ/dV; require a sign-resolved bracket, rather than relying on a near-zero interpolated value alone.

## Primary sources

- Perdigón-Toro et al., *Understanding the Role of Order in Y-Series Non-Fullerene Solar Cells to Realize High Open-Circuit Voltages*, Advanced Energy Materials (2022), DOI: https://doi.org/10.1002/aenm.202103422; main Table 1, Figures 1d, 3, 5; SI Note 2, Note 4, Table S3, Figures S12, S15, S16. The provided original main/SI PDFs were used alongside corrected extraction checks.
- Bässler, *Charge Transport in Disordered Organic Photoconductors: A Monte Carlo Simulation Study*, physica status solidi (b) 175 (1993) 15–56, https://doi.org/10.1002/pssb.2221750102
- Pasveer et al., *Unified Description of Charge-Carrier Mobilities in Disordered Semiconducting Polymers*, Physical Review Letters 94 (2005) 206601, https://doi.org/10.1103/PhysRevLett.94.206601

## Final selected spatial and root checks

With final relative-QF coordinates, local continuity L1 auditing and explicit low-temperature mobility-extrapolation opt-in:

- All 12 selected N641 illuminated states at 100/200/300 K pass current, Poisson charge and chemical-work checks
- Jsc changes on refining N321→641 are 0.2197%, 0.02165%, 0.002586% at 100/200/300 K, respectively
- One further N1281 short-circuit solve at 100 K gives Jsc=−3.4960×10⁻⁹ A cm⁻², only 0.0568% from N641, consistent with approximately second-order spatial convergence
- Independently solved, sign-resolved N641 Voc brackets are 1.1046326–1.1046387 V (100 K), 0.9294360–0.9294421 V (200 K), and 0.7996353–0.7996415 V (300 K)
- Central N321→641 root changes are only 2.99, 1.14 and 0.17 microvolts; no spatial Voc difference is resolved at approximately10-microvolt root precision
- The 100 K bracket required extending the original1.1 V scan; bulk PIA-to-QFLS≈0.9935 V is not the conditional full-device Voc

The extremely small 100 K photocurrent is therefore a resolved consequence of the assumed mobility extrapolation and contact scenario, rather than an unresolved spatial mesh artifact. It must not be reported as an experimentally validated PM6:Y6 prediction.

`FINAL_VALIDATION.json` re-evaluates the selected refined states and root endpoints against final source hashes. `saved_state_audit.json` separately audits main-directory trial outputs: some preliminary dark points fail the tightened final current certificate. These must be re-solved or labeled/excluded by the producer before release; they do not invalidate the independently certified selected illuminated states.
