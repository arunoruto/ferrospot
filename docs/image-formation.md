# Image formation: making the simulated images realistic

Survey of 2026-10-16 for the forward model in `ferrospot.simulate`.
What the model does today, the physics it approximates, and the upgrades
the literature supports, in the order they pay off. Numbers marked
*computed* were calculated for this page (Lenz cross-sections, abTEM
scattering factors, a Monte Carlo of Bragg orientations, abTEM timings on
a 16-core CPU); scripts are in the session scratchpad, not in the
repository. Papers not read in full are listed at the end.

## The current model, and where its numbers come from

Per pixel the electron path crosses $L$ nm of iron (the sum of the
sphere chords) and $t - L$ nm of matrix,

$$
\tau = k_m (t - L) + k_\text{Fe} L, \qquad
I_\text{BF} = e^{-\tau}, \qquad
S_\text{HAADF} = 1 - e^{-\tau},
$$

then Gaussian blur (FWHM 2 px), binning, percentile stretch, white
Gaussian noise and 8-bit quantisation. Only the ratio $k_\text{Fe}/k_m$
matters after the stretch: 4 for HAADF, 3 for bright field.

**Bright field: mass-thickness contrast** {cite:p}`reimer2008tem`.
Electrons scattered elastically beyond the objective aperture
semi-angle $\alpha$ are lost, giving Beer–Lambert in mass thickness with
the screened-Rutherford (Lenz) cross-section,

$$
\frac{I}{I_0} = \exp\Big(-\sum_j n_j \sigma_j(>\alpha) L_j\Big), \qquad
\sigma(>\alpha) = \frac{Z^{4/3}\lambda^2 (1 + E/E_0)^2}{\pi}\,
\frac{1}{1 + \alpha^2/\theta_0^2}, \quad
\theta_0 = \frac{\lambda}{2\pi a_H Z^{-1/3}},
$$

with $n_j$ the atom density, $a_H$ the Bohr radius and $E_0 = 511$ keV.
The Lenz total elastic mean free path reproduces the measured value for
amorphous SiO$_2$ at 200 kV, 178 nm against 172 ± 10 nm (elastic,
measured with a 2 mrad aperture; the inelastic one is 178 ± 6 nm)
{cite:p}`basha2022mfp`, so the ratios are defensible. *Computed* at
200 kV (300 kV changes the ratios by under 5 %):

| aperture $\alpha$ | $k_\text{Fe}$ (nm$^{-1}$) | $k_\text{Fe}/k_m$: SiO$_2$ / MgSiO$_3$ glass / (Mg$_{0.8}$Fe$_{0.2}$)SiO$_3$ |
|---|---|---|
| total | 0.0254 (MFP 39 nm) | 4.5 / 3.1 / 3.0 |
| 5 mrad | 0.0241 | 4.7 / 3.2 / 3.1 |
| 10 mrad | 0.0211 | 5.1 / 3.6 / 3.3 |
| 20 mrad | 0.0141 | 6.1 / 4.3 / 3.9 |

So the current 3 is at the low end; 3 to 4.5 covers silicate glasses at
the usual apertures. In physical units, 75 nm of matrix at 10 mrad gives
$\tau \approx 0.45$ and a 20 nm iron chord adds 0.33: not small, so keep
the exponential rather than linearising. Unfiltered inelastic scattering
passes the aperture and only adds chromatic blur
{cite:p}`iakoubovskii2008`.

**HAADF: incoherent imaging** {cite:p}`pennycook1991,hartel1996`. The
image is the probe intensity convolved with an object function,
$I = |P|^2 \otimes O$, $O \approx \sum_j n_j \sigma^H_j L_j$, with
$\sigma^H$ the cross-section into the detector's angular range,
empirically $\propto Z^n$. Treacy {cite:p}`treacy2011` tabulates the
elastic exponent against the inner collection vector $q_0$ (at 200 kV,
$q_0$ = 1, 2, 4 Å$^{-1}$ is about 25, 50, 100 mrad): for Fe 1.77, 1.88,
1.95 and for Mg and Si within 0.02 of that, so for the matrix-to-iron
contrast the exponent is 1.85 to 1.95 at the usual inner angles, and the
lower measured values (1.64 for single atoms {cite:p}`krivanek2010`, 1.3
to 1.9 {cite:p}`yamashita2018`) come from the inelastic contribution,
which Treacy corrects for with a simple additive term. *Computed*
from Lobato scattering factors at 200 kV (no Debye–Waller or thermal
diffuse scattering, so upper bounds): effective $n$ between O and Fe
1.79 at 30 to 100 mrad, 1.91 at 60 to 200 mrad; $k_\text{Fe}/k_m$ at 60
to 200 mrad 7.2 for SiO$_2$, 5.1 for MgSiO$_3$, 4.4 for an Fe-bearing
glass; with $n$ = 1.6 to 1.8 the ranges are 5.7 to 6.8, 4.0 to 4.7 and
3.6 to 4.2. The current 4 is at the low end. Absolute values
$k_\text{Fe} \approx 1/430$ nm and $k_m \approx 1/2200$ nm: HAADF is
linear in mass-thickness at these thicknesses, and the integrated
intensity of a particle grows monotonically and close to linearly with
its atom count up to 6500 atoms, $R$ = 0.997 for size-selected Au
clusters {cite:p}`young2008`. Validation should compare
particle-integrated intensities, which do not depend on probe size,
defocus or source size {cite:p}`e2013pics`, not peak intensities.
Non-linear damping sets in for long heavy chords
{cite:p}`vandenbroek2012`; Zhong et al. {cite:p}`zhong2018` use exactly
the form $p = I_0[1 - \exp(-\sum_e \mu_e t_e)] + p_b$ and note that
amorphous materials satisfy the projection assumption well in
multislice. Aveyard et al. {cite:p}`aveyard2017` show this by multislice
(200 kV, 10 mrad probe, 100 to 250 mrad detector): amorphous materials
follow a Beer–Lambert-like, nearly linear projection relation, crystals
on a zone axis a severely non-linear one from channelling, and the
linearity improves with higher voltage, a larger inner angle (150 better
than 100 mrad) and a larger convergence angle, and for crystals tilted
10° off axis. For ferrospot this means randomly oriented bcc Fe spheres
have an orientation-dependent HAADF brightness, so a per-particle
visibility factor is justified in HAADF too, not only in bright field
(`BRAGG["haadf"]` is 0 today). No paper was found that compares a chord-length model with
multislice for nanoparticles embedded in glass, so the abTEM validation
below would be new.

**Probe and beam broadening.** The probe is
$\psi_p = \mathcal F^{-1}[A_\alpha(k) e^{-i\chi(k)}]$ with $\chi = \pi
\lambda \Delta f k^2 + \frac{\pi}{2} C_s \lambda^3 k^4$, convolved with a
Gaussian source; at 200 kV and 20 mrad the probe is about 0.15 nm, far
below any pixel here. What is not negligible is beam broadening through
the section: Reed's estimate $b = 7.21 \times 10^5 (Z / E_0)\sqrt{\rho /
A}\, t^{3/2}$ cm gives 1.6 nm at 50 nm and 4.6 nm at 100 nm depth for a
silicate at 200 kV (90 %-current widths, so an overestimate of the blur)
{cite:p}`dejonge2010`. Particles deep in the section are blurrier than
those at the entrance face, which the model ignores.

**Phase contrast in bright field** {cite:p}`kirkland2020`. The exit
wave is $\psi = e^{i\sigma V_p - \mu/2}$ with the projected potential
$V_p$ and interaction constant $\sigma$ = 0.00729 rad V$^{-1}$ nm$^{-1}$
at 200 kV; for a weak phase object the image spectrum is $\delta(k) +
2\sigma \tilde V_p(k) \sin\chi(k) E_s E_t$ with the temporal and spatial
envelopes $E_t = e^{-\frac12 \pi^2 \lambda^2 \Delta^2 k^4}$, $E_s =
e^{-(\pi\alpha_i/\lambda)^2 (\Delta f \lambda k + C_s \lambda^3 k^3)^2}$.
Mean inner potentials (*computed* from neutral atoms, bonding lowers them
by 5 to 10 %): Fe 29 V, SiO$_2$ glass 10.4 V (measured 10.1 ± 0.6 V
{cite:p}`wang1997mip`), MgSiO$_3$ 15.7 V. A 10 nm iron sphere shifts the
phase by about 0.9 rad: a strong phase object, so if phase contrast is
modelled, propagate $e^{i\phi}$ exactly. Fresnel fringes have width
$\sqrt{\lambda \Delta f}$, 1.1 nm at 0.5 µm and 2.2 nm at 2 µm defocus;
underfocus puts a bright fringe outside a dense particle, and the
delocalisation biases apparent sizes by several percent
{cite:p}`yao2012`. The speckle of an amorphous matrix is a Gaussian
random field with power spectrum $\propto t \sum_j n_j |f^e_j(k)|^2$
filtered by $\sin\chi$, with grain size set by the first CTF maximum
$k_1 \approx (2\lambda\Delta f)^{-1/2}$. That is the "phase-contrast
mottle" the calibration cannot reproduce for the Pieters & Noble panels.

**Diffraction contrast** {cite:p}`williams2009tem`. In two-beam
dynamical theory a crystal of thickness $t$ diffracts

$$
I_g = \Big(\frac{\pi t}{\xi_g}\Big)^2 \frac{\sin^2(\pi t s_\text{eff})}{(\pi t s_\text{eff})^2},
\qquad s_\text{eff} = \sqrt{s^2 + \xi_g^{-2}}, \qquad
I_\text{BF} = e^{-\tau}\Big(1 - \sum_g I_g\Big),
$$

with extinction distances for bcc Fe at 200 kV (*computed*) of 36 nm
(110), 52 nm (200), 67 nm (211) and 82 nm (220), and $2\theta_B(110)$
= 12.4 mrad, so Bragg contrast appears only when the aperture is smaller
than that. A Monte Carlo over random orientations with $t = D$
(*computed*) gives the share of particles losing more than 5 / 20 / 50 %
of the bright-field beam:

| $D$ (nm) | > 5 % | > 20 % | > 50 % |
|---|---|---|---|
| 2 | 2 % | 0 | 0 |
| 5 | 23 % | 1 % | 0 |
| 10 | 20 % | 9 % | 3 % |
| 20 | 17 % | 8 % | 5 % |
| 40 | 17 % | 6 % | 2 % |

Roughly 5 to 10 % of 10 to 30 nm particles are black and 10 to 15 %
grey, particles under 3 nm show none, and about 80 % show no significant
diffraction contrast in the one literature count found. The present
`BRAGG` defaults (25 % strong at 4×, 25 % mid at 2×) are too many and
too strong; because $t = L(x)$ varies across a sphere, diffracting
particles also show thickness fringes and partial, ring-like darkening.

## Noise and detectors

The signal is Poisson in the dose, $N \sim \text{Poisson}(D a^2 I)$ with
$D$ electrons per area and $a$ the pixel; HAADF collects 1 to 10 % of the
probe, so it is Poisson-limited. A scintillator camera spreads each
electron, so the detector MTF applies *after* the Poisson draw and the
noise becomes correlated with spectrum $\propto |\text{MTF}|^2$; read
noise and gain come last {cite:p}`ruskin2013`. Scintillator CCDs have
DQE(0) of 0.3 to 0.5 falling to about 0.1 at Nyquist
{cite:p}`mcmullan2009,mcmullan2014`. Networks trained on simulations
only transfer to low-dose experiments when the camera MTF is modelled;
the exact MTF is not needed, a range of Lorentzian-like profiles
(parametrised by its half-maximum frequency and curvature) works, and
the noise matters most below 200 e$^-$/Å$^2$ with detection still
possible at 20 to 30 e$^-$/Å$^2$ {cite:p}`lethlarsen2023`. STEM adds scan noise: per-line jitter, slow
drift and flyback distortion {cite:p}`jones2013,ophus2016drift,
sang2014`. Published figures add resampling, contrast stretch and
gamma, unsharp masking, 8-bit quantisation and JPEG blocking; no paper
treats these, so they are augmentations.

Order of operations:

| step | TEM (BF, DF) | STEM (HAADF) |
|---|---|---|
| 1 | object: $\tau$, or exit wave | object |
| 2 | aperture or CTF, intensity | probe ⊗ object ⊗ source |
| 3 | Poisson | Poisson |
| 4 | detector MTF, read noise | scan jitter, drift, amplifier noise |
| 5 | resample, stretch, quantise, JPEG | same |

## Multislice: the reference the cheap model is checked against

Multislice {cite:p}`cowley1957` propagates $\psi_{n+1} = p_{\Delta z}
\otimes [t_n \psi_n]$, $t_n = e^{i\sigma v_n}$, $\tilde p = e^{-i\pi
\lambda \Delta z k^2}$, with frozen phonons ($u_\text{rms} = \sqrt{B /
8\pi^2}$, 0.07 Å for Fe). Packages: abTEM {cite:p}`madsen2021abtem`
(Python, NumPy/CuPy, Dask, ASE; `pip install abtem` works CPU-only,
version 1.0.10 verified), PRISM and Prismatic
{cite:p}`ophus2017prism,pryor2017,rangeldacosta2021`, Dr. Probe
{cite:p}`barthel2018`, MULTEM {cite:p}`lobato2015,lobato2016`. abTEM has
no glass builder: place atoms at random at the right density (96
atoms/nm$^3$ for MgSiO$_3$ at 3.2 g/cm$^3$) with a minimum-distance
rejection, cut a sphere and insert a randomly rotated bcc Fe cluster
from ASE. The calls are `Potential(atoms, sampling, slice_thickness,
parametrization="lobato")`, `FrozenPhonons`,
`PlaneWave(energy).multislice(pot).apply_ctf(defocus, Cs).intensity()`
for BF, `Probe(energy, semiangle_cutoff).scan(pot, GridScan, AnnularDetector(60, 200))`
for HAADF, `SMatrix(interpolation=f)` for PRISM, and `Images` carries
`poisson_noise(dose_per_area)`, `gaussian_filter`, `scan_noise`.

Measured on a 16-core CPU at 200 kV, 2 Å slices, a 10 × 10 × 20 nm cell
of 190 000 atoms with one 3 nm-radius Fe sphere:

| run | grid | potential | multislice |
|---|---|---|---|
| plane wave, 0.1 Å sampling | 1000² | 14 s | 7 s |
| plane wave, 0.05 Å sampling | 2000² | 51 s | 30 s |
| HAADF, conventional, 0.05 Å | 2000² | | 8 s per probe position |

PRISM with $f = 4$ did not finish in 20 min on CPU and took 18 GB.
Extrapolated: a 30 × 30 × 100 nm bright-field cell is about 15 min per
frozen-phonon configuration at 0.1 Å (potential 72 GB unless streamed);
HAADF of 10 × 10 × 100 nm at 0.5 nm/px is 400 positions × 40 s, 4 to 5 h
per configuration on CPU. So: bright field on CPU in minutes, HAADF on
small cells with GPU PRISM (10 to 50× faster), never the calibration
grid of 200 images per panel.

**Validation protocol.** Frozen phonons, random-position glass,
randomly rotated bcc spheres in 10 × 10 × (50 to 100) nm cells. Compare
with the cheap model: particle-integrated intensity against diameter
and depth {cite:p}`e2013pics`, the background power spectrum, and edge
profiles across particles. Fit from these the free parameters of
upgrades 1, 2 and 6 below.

## Synthetic data for detectors

The HRTEM framework of Rangel DaCosta et al. {cite:p}`rangeldacosta2024`
(Prismatic, 8 frozen phonons, aberrations, focal spread, plasmon
losses, Poisson at 200 to 600 e$^-$/Å$^2$, no MTF) validates only
indirectly, by detector F1 on experimental benchmarks, and finds that
the largest gains come from a wide spread of noise levels and substrate
thicknesses, not from more structures. Madsen et al.
{cite:p}`madsen2018` trained purely on multislice and transferred to
experiment; Vincent et al. {cite:p}`vincent2021` did the same for
denoising. None compares pixel statistics of simulated and real images
directly, which is the check ferrospot should add: integrated intensity
against size, background power spectrum, edge profiles, and the
pixel-scale noise residual already matched in `calibrate.match_noise`.

## Upgrades, in order of realism per unit cost

1. **Physical coefficients instead of fixed ratios.** $k^\text{BF} =
   \sum_j n_j \sigma_j(>\alpha)$ from Lenz, $k^H \propto \sum_j n_j
   Z_j^n$ with $n \sim U(1.6, 1.9)$; randomise $k_\text{Fe}/k_m$ over
   [3, 4.5] (BF) and [3.6, 5.5] (HAADF) across the grid as a nuisance,
   and keep $\tau$ in physical units (matrix $\tau_\text{BF}$ 0.3 to
   0.6) so the exponential saturation is right.
2. **Depth-dependent blur.** Each sphere at depth $z$ gets
   $\text{PSF}_z = e^{-z/\ell}|P|^2 + (1 - e^{-z/\ell})\, G(r; w \approx z
   \bar\theta)$ with $\ell$ = 110 to 170 nm and $\bar\theta \approx 20$
   mrad: a 1 to 2 nm halo on deep particles, which matters at pixels
   below 1 nm. Cheap: bin spheres by depth, blur each bin's projection
   separately.
3. **Two-beam Bragg term along the chord** replacing the random
   visibility factor: $I_\text{BF} = e^{-\tau}[1 - \sum_g I_g(L(x),
   s_g)]$ over {110, 200, 211} with a random orientation per sphere and
   the extinction distances above, applied only when $\alpha <
   2\theta_B$. Gives the right black and grey fractions, the size trend
   and internal fringes, with no FFT.
4. **Poisson, then MTF, then read noise** instead of white Gaussian
   noise: $N = \text{Poisson}(D a^2 I)$, $N' = \mathcal F^{-1}[\text{MTF}
   \cdot \mathcal F N] + \mathcal N(0, \sigma_r^2)$, with an empirical
   MTF such as $\text{sinc}(k / 2k_N)[a e^{-k^2/2b^2} + (1 - a)/(1 +
   (k/c)^2)]$, parameters and dose randomised. The noise residual then
   has the correlated spectrum real figures have.
5. **Exit wave and CTF for bright field at ≤ 0.5 nm/px** (two FFTs):
   $\psi = e^{i\sigma[V_m(t - L) + V_\text{Fe} L + \delta V_\text{sp}] -
   \tau/2}$, $I = |\mathcal F^{-1}[\tilde\psi A_\alpha e^{-i\chi} E_s
   E_t]|^2$ with $V_\text{Fe}$ = 25 to 29 V, $V_m$ = 10 to 15 V, $\Delta
   f \sim U(-2, +1)$ µm and $\delta V_\text{sp}$ the speckle field above.
   This is what the Pieters & Noble panels need; at coarse pixels it
   reduces to upgrade 1.
6. **HAADF saturation**, $S = \frac{\sum_j k^H_j L_j}{\sum_j k^\text{tot}_j
   L_j}[1 - e^{-\sum_j k^\text{tot}_j L_j}]$ with total elastic
   coefficients (Fe 1/28 to 1/39 nm, matrix 1/110 to 1/170 nm), a
   heuristic consistent with Zhong 2018, to be calibrated against abTEM.
7. **STEM scan artefacts**: row shifts $\sim \mathcal N(0, \sigma_j^2)$
   with $\sigma_j$ = 0.1 to 0.5 px, affine drift, flyback.
8. **Published-figure augmentation**: anti-aliased resampling with a
   random kernel, percentile clip 0.1 to 2 % and gamma 0.8 to 1.25,
   optional unsharp mask, 8-bit, JPEG quality 70 to 95.
9. **Probe and partial coherence** only below about 0.3 nm/px.
10. **The abTEM validation** above, as an optional `multislice` extra
    and a script, not part of the calibration loop.

## Not read in full

Read in full (local copies): Treacy 2011, Young 2008, Aveyard 2017,
Leth Larsen 2023, Basha 2022. Abstract only: Van den Broek 2012
(10.1016/j.ultramic.2012.03.005), Yildirim & Cole 2021
(10.1021/acs.jcim.0c01455), Yao 2012 (10.1093/jmicro/dfs065,
unverified), Wang 1997 (10.1063/1.118556, unverified). Not retrieved:
Iakoubovskii 2008 tables (10.1002/jemt.20597).
Cited from memory with DOIs to check: the three textbooks, Pennycook &
Jesson, Hartel, Krivanek, both LeBeau papers, Rosenauer, Van Aert, E et
al., Cowley & Moodie, the multislice package papers, Ruskin, both
McMullan papers, Jones & Nellist, Ophus 2016, Sang & LeBeau, Madsen 2018,
Vincent 2021.
