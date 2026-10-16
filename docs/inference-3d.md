# From a projection to 3D: what the literature offers

Survey of 2026-10-16 for the inverse problem ferrospot solves by
simulation: a FIB section of thickness $t$ holds spheres of diameter $D$
at number density $n_V$; the image shows their projection. Equations
marked *derived* were worked out for this page from standard stochastic
geometry and checked numerically; they are not quoted from a paper.
Papers that could not be read in full are listed at the end.

## Counts: projected area to volume

**Cahn–Nutting** {cite:p}`cahn1959`: a foil shows every sphere that
touches it, so

$$
N_A = n_V\,(t + \bar D),
$$

with $\bar D$ the *true* mean diameter. ferrospot uses the observed
half-maximum diameter here, which is smaller (caps cut by the faces,
the $\sqrt 3/2$ outline), so $n_V$ comes out slightly high. Hilliard
{cite:p}`hilliard1962` adds a visibility cutoff: a cap thinner than
$h_c$ is not seen, so $N_A = n_V (t + D - 2 h_c)$ (*derived* form for a
thickness-contrast cutoff). The Holmes effect {cite:p}`underwood1970`:
opaque particles inflate the projected area fraction, $A_A = V_V +
\tfrac14 S_V t$, for spheres $A_A = V_V (1 + 3t / 2D)$.

**Observed sizes in a thick foil** (*derived*, same kernel as Goldsmith
{cite:p}`goldsmith1967`): a fraction $t/(t + D)$ of the spheres lie wholly
inside and show $D$; the rest are caps with chord-diameter density
$f(d) = d / (D\sqrt{D^2 - d^2})$, the Wicksell kernel
{cite:p}`wicksell1925`, so

$$
\mathbb E[d_\text{obs}] = D\,\frac{t + \pi D / 4}{t + D}, \qquad
N_A(d) = t\, n_V(d) + d \int_d^\infty \frac{n_V(D)}{\sqrt{D^2 - d^2}}\, dD .
$$

$t = 0$ is Wicksell's planar-section problem (Saltykov unfolding); for
$t \gg D$ the first term dominates and the foil is nearly unbiased for
size. For $t = 100$ nm and $D = 3$ nm the cap bias is below 1 %, so size
unfolding is not worth it for the SnpFe0 rims; it matters for the 30 nm
LnpFe0.

**Hidden overlaps** {cite:p}`miles1976`: for overlapping (Boolean)
grains in a slab the projected area fraction is

$$
A_A = 1 - \exp\!\left[-n_V\left(\bar V + t \bar A\right)\right], \qquad
\bar V = \tfrac43 \pi\, \mathbb E[R^3],\; \bar A = \pi\, \mathbb E[R^2],
$$

which inverts to a *count-free* density estimate $n_V = -\ln(1 - A_A) /
(\bar V + t \bar A)$ that already accounts for particles hidden behind
others. In the dilute limit it is the Cahn–Nutting formula. The expected
number of projected overlaps per particle for hard spheres of exclusion
$D$ (*derived*) is

$$
N_\text{ov} = n_V \left[\pi D^2 T - \tfrac43 \pi D^3 + \frac{\pi D^4}{2T}\right],
\qquad T = t + D,
$$

the Poisson term minus the hard-core reduction ($= 8\phi$); the share of
particles with at least one projected overlap is about $1 -
e^{-N_\text{ov}}$. That is an analytic check on the overlap fraction the
pipeline measures.

**Edge weights** {cite:p}`miles1974`: a particle wholly inside an $a
\times b$ window with bounding box $w \times h$ gets weight $ab / [(a -
w)(b - h)]$; large particles touch the border more often, and this makes
size and count distributions unbiased.

## Nearest neighbours

**2D Poisson** {cite:p}`ClarkEvans1954`: NN density $p(r) = 2\pi\lambda r
e^{-\pi \lambda r^2}$, mean $1 / (2\sqrt\lambda)$, variance $(4 -
\pi) / (4\pi\lambda)$, standard error of the ratio $0.26136 /
\sqrt{N\lambda}$. Donnelly's edge correction for a rectangle of area $A$
and perimeter $P$ {cite:p}`donnelly1978`:

$$
\mathbb E[\bar r] \approx 0.5\sqrt{A/N} + \left(0.0514 + \frac{0.041}{\sqrt N}\right)\frac{P}{N}, \qquad
\operatorname{Var} \approx 0.0703\,\frac{A}{N^2} + 0.037\, P \sqrt{A / N^5}.
$$

The border correction ferrospot uses (keep a particle if its NN is
closer than the edge) is biased by $-3$ % ({doc}`optimisation`); the
Hanisch estimator {cite:p}`hanisch1984` and, in 3D, Reed et al.
{cite:p}`reed1997` are the unbiased alternatives for the full NN
distribution $G(r)$.

**3D Poisson** {cite:p}`hertz1909,Chandrasekhar1943`: $p(r) = 4\pi n r^2
e^{-\frac43 \pi n r^3}$, $\mathbb E[r^k] = \Gamma(1 + k/3)(4\pi n /
3)^{-k/3}$, so mean $0.554\, n^{-1/3}$ and standard deviation $0.201\,
n^{-1/3}$. Projecting a 3D Poisson process gives a 2D Poisson process
of intensity $\lambda = n_V (t + D)$, so $\mathbb E[r_{3D}] / \mathbb
E[r_{2D}] = 1.108\, n_V^{1/6} \sqrt{t + D}$ (*derived*). That is the
whole content of the Poisson estimate in `spacing_stats`.

**Hard spheres** {cite:p}`torquato1990,torquato1995pre`: for equilibrium
hard spheres of diameter $\sigma$ at volume fraction $\phi$, with $x =
r/\sigma \ge 1$,

$$
E_V(r) = \exp\{-\phi[8a_0(x^3 - 1) + 12a_1(x^2 - 1) + 24a_2(x - 1)]\},
$$
$$
a_0 = \frac{1 + \phi + \phi^2 - \phi^3}{(1 - \phi)^3},\quad
a_1 = \frac{\phi(3\phi^2 - 4\phi - 3)}{2(1 - \phi)^3},\quad
a_2 = \frac{\phi^2(2 - \phi)}{2(1 - \phi)^3},
$$

is the probability that the nearest centre is farther than $r$; the NN
density is $-dE_V/dr$ and the mean $\ell = \sigma + \int_\sigma^\infty
E_V\, dr$. Against the Poisson value (computed for this page):

| $\phi$ | hard-sphere $\ell / D$ | Poisson $0.554 n^{-1/3} / D$ | ratio |
|---|---|---|---|
| 0.01 | 2.16 | 2.07 | 1.04 |
| 0.05 | 1.41 | 1.21 | 1.16 |
| 0.10 | 1.22 | 0.96 | 1.27 |
| 0.20 | 1.10 | 0.76 | 1.44 |

At 5 % fill the Poisson surface gap is $0.21 D$ against $0.41 D$ for
hard spheres: the Poisson readout halves the gap, which is why
`spacing_stats` calls it a rough guide and the simulation replaces it.
Random sequential addition (what polydisc does) differs little from the
equilibrium ensemble below $\phi \approx 0.2$ {cite:p}`torquato2006rsa`;
it jams at 0.384.

**Projected NN of hard spheres**: no closed form found. Li et al.
{cite:p}`li2014disp` project hard-core processes through 50 to 100 nm
sections and show the projection tends to 2D Poisson only as $t \to
\infty$; at realistic thickness, Poisson-based (CSR) metrics mislead, and
Ripley's $K$ was the metric most robust to a wrong thickness and most
powerful against clustering.

## Pair correlation of a projected slab

Centres occupy an effective slab $T = t + D$ at density $n_V$. The
projected pair correlation (*derived*) is

$$
g_2(s) = \frac{2}{T} \int_0^T \left(1 - \frac{z}{T}\right)
g_3\!\left(\sqrt{s^2 + z^2}\right) dz
\;\approx\; 1 + \frac{1}{T}\, \mathcal A[g_3 - 1](s), \qquad
\mathcal A[h](s) = 2\int_s^\infty \frac{h(r)\, r\, dr}{\sqrt{r^2 - s^2}},
$$

an Abel transform when $T$ exceeds the correlation length, with inverse
$h_3(r) = -\frac1\pi \int_r^\infty \frac{d[T h_2]/ds}{\sqrt{s^2 - r^2}}\,
ds$. For a hard core of exclusion $\sigma$ at low density, with $u =
\sqrt{\sigma^2 - s^2}$, $g_2(s) = 1 - \frac{2u}{T}(1 - \frac{u}{2T})$ for
$s < \sigma$: the hole is only $2\sigma / T$ deep ($g_2(0) \approx 0.88$
for $D = 5$ nm, $T = 80$ nm). Hard-core order is diluted by $D / T$ in
projection and needs many particles to see. Clustering is not: a 3D
Thomas process {cite:p}`thomas1949` (parents $\kappa$, Gaussian spread
$s_c$) projects to a 2D Thomas process with parents $\kappa T$ and the
same $s_c$,

$$
g_2(s) - 1 = \frac{1}{\kappa T\, 4\pi s_c^2}\, e^{-s^2 / 4 s_c^2},
$$

so $g_2(s)$ or Besag's $L(s) = \sqrt{K(s)/\pi}$ separate clustered from
hard-sphere arrangements far better than the mean NN distance.

## Ground truth: tomography and atom probe

* Dai et al. 2026 {cite:p}`dai2026`: HAADF electron tomography (±80°
  tilt, Titan Themis at 300 kV) of two FIB tips from the impact-glass
  coating of a Chang'e-5 grain, plus a 70 nm lamella. The npFe0 are
  highly spherical bcc particles; only those above 25 nm in the
  large-particle plane are irregular, with necks from sintering. Per
  tomogram:

  | tomogram | particles | mean $D$ | median $D$ | Fe$^0$ volume fraction |
  |---|---|---|---|---|
  | Tomo-1 (six layers) | 1506 | 3.4 nm | 2.9 nm | 2.9 % |
  | Tomo-2 | 43 | 5.9 nm | 5.0 nm | 2.8 % |
  | Tomo-3 | 385 | 8.3 nm | 7.9 nm | 2.0 % (2.8 % outside the particle-free glass) |
  | Tomo-4 | 18 | 12.7 nm | 11.6 nm | 2.8 % |

  Within Tomo-1 the layers L1 and L3 to L6 hold 2 to 3 vol % with mean
  diameters rising from 3.2 nm (L1) to 5.0 nm (L6, nearest the surface),
  and the large-particle plane L2 reaches 30 vol %. No nearest-neighbour
  spacing is reported, but the particle list exists (Fig. 4A), so the
  authors could be asked for it. What the numbers imply: at 2.8 vol %
  and $D$ = 3.4 nm the number density is about $1.4 \times 10^6$
  µm$^{-3}$, and the hard-sphere mean NN centre distance ($\ell / D
  \approx 1.6$ at $\phi$ = 0.03 from the Torquato formula above) is
  about 5.5 nm, a mean surface gap of roughly 1.2 radii. That is denser
  than ferrospot's calibrated estimate for the Shen 2024 SnpFe0 rim (0.7
  vol %, gap 2 to 8 radii) but a different sample, an impact-glass
  coating rather than a solar-wind rim; it is the right order of
  magnitude to test the pipeline on, by projecting the tomography
  particle list through 70 nm and comparing with Slice-1.
* Greer et al. 2026 {cite:p}`greer2026`: TEM and atom-probe tomography on
  Apollo 17 soil; npFe0 mostly 3 to 5 nm; APT cluster counting "fails to
  meaningfully characterize" the population; early npFe0 form a
  sponge-like network of Fe-rich clusters; 2D projection artefacts named
  as a limitation.
* Tomography of supported catalyst particles with 3D NN distances
  {cite:p}`ward2007,prieto2013,friedrich2009`: the method transfers.

## npFe0 in the literature

Beyond the papers already at hand (Keller & McKay 1997, Pieters & Noble
2016, Thompson 2016 and 2017, Li 2022, Shen 2024): Noble, Pieters &
Keller 2007 {cite:p}`noble2007` give synthetic-npFe0 size distributions
and Fe wt %, convertible to a volume fraction with $\phi = (w /
\rho_\text{Fe}) / (w / \rho_\text{Fe} + (1 - w) / \rho_m)$; Burgess &
Stroud 2018 {cite:p}`burgess2018` find npFe0 below 3 to 4 nm in the
silicate rims of one Apollo grain, but in its ilmenite rim elongated,
vesicle-associated inclusions aligned with the host lattice, so sphere
models do not transfer to oxide hosts; Dai et al. 2025
{cite:p}`dai2025` describe clustered, partly merged Fe-rich nanodroplets
in impact glass; Gu 2022 {cite:p}`gu2022`, Guo 2022 {cite:p}`guo2022`
and Chaves & Thompson 2022 {cite:p}`chaves2022` give sizes but no
spacings. No paper found reports a number density or NN spacing except
Dai 2026. The optical motivation: when the surface gap is about one
diameter or less, near-field coupling changes absorption relative to the
isolated-sphere assumption, computable with MSTM
{cite:p}`mackowski2011`.

## Sizing and detection

Lindeberg {cite:p}`Lindeberg1998`: the scale-normalised LoG peaks at
$\sigma^\star = r / \sqrt 2$ for a flat disc; for the projected-sphere
chord profile it is $\sigma^\star = 0.602 R$ (numerical, matching
`SIGMA_PER_RADIUS`), so $R = 1.66\, \sigma^\star$. The generalised LoG of
Kong et al. {cite:p}`kong2013` handles elliptical and touching blobs.
The FWHM of a projected sphere is $\sqrt 3 R$ unblurred, but with a
Gaussian PSF of standard deviation $\sigma_p$ the ratio FWHM$/D$ is not
monotonic: 0.87, 0.84, 0.80, 0.86, 1.30 for $\sigma_p / R$ = 0, 0.1, 0.25,
0.5, 1 (computed for this page), so a single correction factor does not
exist for 1 to 3 nm particles and the forward model has to carry the
bias, as the calibration does. Pyrz & Buttrey {cite:p}`pyrz2008` is the
practical reference for pixels per particle, threshold and defocus
sensitivity. Deep-learning segmentation {cite:p}`groschner2021,
yildirim2021,larsen2023npsam` separates touching instances, NP-SAM with an
explicit treatment of partial overlap, but none models projection
overlap of 3D particles physically.

## A third packing family: clustered

Thomas {cite:p}`thomas1949` and Matérn cluster processes
{cite:p}`chiu2013` have closed-form $g(r)$; fractal aggregates
{cite:p}`meakin1983,kolb1983,lin1989` scale as $N = k_f (R_g / a)^{D_f}$
with $D_f \approx 1.8$ (DLCA) to 2.1 (RLCA), and FracVAL
{cite:p}`Moran2019` builds them with tunable $D_f$, $k_f$ and
polydispersity. The stronger evidence in the rims is depth
stratification (Dai 2026) and Fe-rich clusters (Greer 2026) rather than
chains, which argues for inhomogeneous statistics against distance from
the grain surface {cite:p}`baddeley2000` first, and a cluster family
second.

## What to add, in order

1. Read the 3D NN distance off the simulated packing (already done in
   `simulate.truth`) and retire the Poisson readout, or replace it with
   the hard-sphere $\ell$ above.
2. Cahn–Nutting with the true $D$ and a cap cutoff.
3. The projected pair correlation $g_2(s)$ and $L(s)$, edge-corrected,
   as matching statistics beside the NN distribution.
4. Miles' area-fraction density as a count-free cross-check, and the
   analytic overlap count as a sanity check.
5. The full edge-corrected $G(r)$ (Hanisch) and Clark-Evans with
   Donnelly's correction and standard error.
6. Depth-resolved statistics, then a Thomas or FracVAL cluster family.

## Not read in full

Read in full (local copies): Dai 2026, Burgess & Stroud 2018. Abstract
only: Chaves & Thompson 2022 (10.1186/s40623-022-01683-6), NP-SAM
(10.26434/chemrxiv-2023-k73qz). Cited from memory with the DOI checked:
Noble 2007, Gu 2022, Guo 2022, Miles 1976, Goldsmith 1967, Torquato 1990,
Ward 2007, Prieto 2013, Pyrz 2008, Groschner 2021, Yildirim 2021. No DOI:
Cahn & Nutting 1959, Hilliard 1962, Underwood 1970, Donnelly 1978, Miles
1974.
