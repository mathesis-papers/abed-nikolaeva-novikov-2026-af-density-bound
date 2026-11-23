# Sharp accretive modulus and density bounds

Haobo Ma (ChronoAI Pte Ltd, Singapore), then Wenlin Zhang (National University
of Singapore; ChronoAI Pte Ltd, Singapore). 6 November 2026. CC BY 4.0.

## Question and results

The note answers the density-independent boundedness question in
Abed–Nikolaeva–Novikov (2026), arXiv:2610.34388v1, Section 11.2(5),
and determines the best dimension coefficient discussed in Remark 5.45.
It applies Hayashi's (2001) pinching inequality and the projection-count
form in Ogawa–Hayashi (2002). Classical triangular models are credited
through Arlinskii–Tsekanovskii (2006), and the dimension-two proof and
fixed-density attainment are credited to Abed–Nikolaeva–Novikov (2026).

For every accretive complex matrix X, Re X <= m(X)|X| <= rank(X)|X|,
where m(X) counts distinct positive singular values. For a positive
definite n-dimensional density of condition number kappa,

    1 <= a(rho) <= n*kappa/(kappa+n-1).

The last expression is below n when n>1. The best dimension coefficient
is n; the supremum over positive definite densities is n, also with
rank-one Hermitian parts. The note includes fixed-rank sharpness,
operator nonattainment at rank greater than one, the positive-imaginary-part
variant and block monotonicity. The coverage table distinguishes these
conclusions from the full density-dependent formula, the II_1 questions,
the numerical data and exact density-multiplicity observations.

The earlier bound a(rho) <= (kappa+1)/(2*sqrt(kappa)) in
Abed–Nikolaeva–Novikov, Proposition 5.44(a), remains applicable; take
the minimum of it and the dimension-dependent upper bound above.
Winter (2025), arXiv:2510.19111v1, Section I, Lemma 1, gives weighted
generalisations of the same classical pinching inequality.

## Literature scope

Moradi–Furuichi–Sababheh, “Norm bounds and factorization results for
accretive matrices,” The Journal of Analysis 34 (2026), 2641–2660,
https://doi.org/10.1007/s41478-026-01123-8, is related prior work.
Full-text access was unavailable; comparison is limited to its abstract
and bibliography. Possible overlap concerns the dimension or distinct-
singular-value Loewner-order bounds, weighted density bound, or sharp
dimension coefficient, including equivalent factorisation consequences.
The sector-angle terminology in the abstract does not rule this out.
The note gives explicit proofs with classical attribution, without a
priority assertion relative to that article.

## Build the PDF

Use Tectonic 0.17.0 with the standard packages listed in `note.tex`:

```sh
tectonic note.tex
```

Alternatively, use a standard LaTeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded; no BibTeX input or external figure is required.

## Exact finite checks

Python 3.9 or newer; standard library only:

```sh
python3 verify.py > checks-new.json
cmp checks-new.json checks.json
```

The deterministic output covers 11 rational matrix examples in dimensions
1–6 (including five singular examples), six examples with repeated singular
values in dimensions 1–6, and eight rational sharpness examples in dimensions
1–8. The sharpness quotients lie within 1/100 of their dimensions. Exact
arithmetic over Q(i) checks the triangular diagonalizers, their inverses,
the positive metric factorization and T*GT=G. Spectral enclosures in the
first group and exact condition numbers in the repeated-singular-value
group test the corresponding quotient bounds. These are finite examples;
the general conclusions follow from the written proofs.
Long integers are displayed as exact sums of powers of ten.

## Lean declarations and written proofs

Pinned source:

https://github.com/the-omega-institute/trureturing/blob/d9bb8a406783c5222849b798f3c2ab738f653d75/D5/S3/Quantum/Matrix/AbedNikolaevaNovikovAFDensityBound.lean

In namespace `D5.S3.Quantum.Matrix.AbedNikolaevaNovikovAFDensityBound`,
`accretive_le_card_abs` proves A <= n|A+iB| for positive semidefinite A
and Hermitian B, including singular matrices. `result` proves nonemptiness,
the least-upper-bound property and the density-independent witness C=n
for the exact quotient set. The singular-value and rank refinements,
condition-number bound, strictness, attainment, positive-B deduction,
block monotonicity and sharpness have written proofs in `note.tex`.
`checks.json` contains the finite examples specified above.

The imported trace definition is `bornProbability rho P = Matrix.trace (rho * P)`:

https://github.com/the-omega-institute/trureturing/blob/d9bb8a406783c5222849b798f3c2ab738f653d75/D5/S3/Quantum/FiniteDimensional.lean

## Source archive

`note.tex` is the sole main source. `metadata.json` gives the creators,
license, statement identity and pinned Lean link. `verify.py` and
`checks.json` contain finite identities and expected output.
`LICENSE.txt` gives the license notice. `SHA256SUMS` identifies the source
members. The readable `note.pdf` accompanies the source archive separately.
