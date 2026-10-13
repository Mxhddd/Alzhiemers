# Network audit (2026-10-13, about 20:50 UTC)

The researcher set cloud network access to "Full". This audit was re-run
from scratch after that change. Earlier results were not assumed.

## Method

- **One plain request per URL.** `curl -sS -L`, default curl user agent,
  40 s limit, first 64 kB only. No retries, no browser headers, no
  workarounds for bot checks, no mirrors or caches for pages that refuse
  automated access. Script: `scratchpad/netaudit/audit.py`, run 2026-10-13;
  its raw JSON output is summarised here.
- **Proxy view.** `curl "$HTTPS_PROXY/__agentproxy/status"` (README:
  `/root/.ccr/README.md`) before and after the probes:
  - `selective: false`, `toolScoped: false`;
  - `recentRelayFailures: []`.
  - `pypi.org` and `files.pythonhosted.org` are in `noProxy`, so they are
    reached directly (CONNECT code 0 below).
- **Reading the columns.**
  - "CONNECT" is the proxy's answer to the tunnel request. 200 means the
    environment allowlist let the request through.
  - "HTTP" is the site's own answer.
  - A blocked host would show CONNECT 403/407, or DNS/TLS errors. None did.
- **No email address was sent anywhere.**

## Result: nothing is blocked by the environment allowlist

Every host was reached: CONNECT 200 (or direct), no DNS errors, TLS valid
with the proxy CA. Every remaining block is the site's own bot protection,
a token requirement, or the session's GitHub repository scoping.

### (a) Blocked by the environment

| Host | What happened | Researcher fix? | Open alternative used |
|---|---|---|---|
| github.com, api.github.com, codeload.github.com, for repositories not attached to this session (tested: `MRCIEU/TwoSampleMR`) | HTTP 403 from the session's GitHub gateway: "GitHub access to this repository is not enabled for this session. Use add_repo" | Not needed | `mrcieu.r-universe.dev` (HTTP 200) serves TwoSampleMR and its dependencies as R packages. `raw.githubusercontent.com` also answers (206). The project's own repository pushes fine. |

### (b) Blocked by the site's bot protection (not bypassed)

All five returned HTTP 403 with Cloudflare's "Just a moment..." page:
`server: cloudflare`, `cf-mitigated: challenge` and the challenge-platform
script.

| Host | URL tested | Open alternative |
|---|---|---|
| academic.oup.com | Braak 2011 article page and PDF | The researcher supplied the PDF. Saved as `data/raw/papers/braak2011.pdf` (git-ignored; sha256 `9bc53256…f2bd2`). |
| jamanetwork.com | home page | Europe PMC REST full text (open-access articles), NCBI PMC, PubMed metadata via E-utilities |
| journals.lww.com | home page | Same as above |
| n.neurology.org | home page | Same as above |
| europepmc.org (website only) | home page | The **Europe PMC REST API** at `www.ebi.ac.uk/europepmc/webservices/rest/` works (search and `fullTextXML`, HTTP 200) |

### (c) Needs an account or token

| Service | What happened | What is needed, and where | Open alternative used |
|---|---|---|---|
| OpenGWAS API (`api.opengwas.io`) and its file downloads (`gwasinfo/files`) | `/api/status` 200, but `/api/user` and `/api/gwasinfo?id=…` return 401: "From 1st May 2024 you must provide a token (JWT)". `gwas.mrcieu.ac.uk/files/…` now redirects to `opengwas.io` and returns 404. File downloads go through the token-gated API | A personal OpenGWAS JWT. Sign in at https://api.opengwas.io/ (account login; this links an email to the token, so only the researcher can do it), open the profile page, generate a token, and store it as an environment secret named `OPENGWAS_JWT` | GWAS Catalog harmonised summary statistics on `ftp.ebi.ac.uk` (200); FinnGen R10 public summary statistics on `storage.googleapis.com` (manifest 200); eQTLGen on `molgenis26.gcc.rug.nl` (200) |
| GWAS Catalog summary-statistics **API** (`www.ebi.ac.uk/gwas/summary-statistics/api/`) | HTTP 410: "This API has been deprecated" | Not a permissions issue; the service was retired | GWAS Catalog FTP over HTTPS, plus REST API v2 (`/gwas/rest/api/v2/…`, 200) |
| CORE API (tested earlier this session, not re-tested) | Needs an API key, issued after registering with an email address | Only the researcher can register. Not needed now: the Braak PDF was supplied | Europe PMC REST, NCBI PMC |

### Reachable (HTTP 200/206; 206 = byte-range answer to the 64 kB limit)

| Host | Use | HTTP | Note |
|---|---|---|---|
| cran.r-project.org, cloud.r-project.org | R packages | 206 | |
| bioconductor.org | Bioconductor | 206 | Behind Cloudflare, but not challenged |
| www.ebi.ac.uk | GWAS Catalog REST v2; Europe PMC REST | 200 | |
| ftp.ebi.ac.uk | GWAS Catalog summary statistics (FTP over HTTPS) | 200 | |
| rest.ensembl.org | gene coordinates | 200 | Slow (14 s) |
| ftp.ensembl.org | archived GTF | 206 | Alternative to REST for gene coordinates |
| eutils.ncbi.nlm.nih.gov | PubMed / PMC metadata | 200 | |
| pmc.ncbi.nlm.nih.gov | PMC article pages | 200 | |
| api.opengwas.io | status endpoint only | 200 | Data endpoints need a token, see (c) |
| raw.githubusercontent.com | raw files | 206 | |
| mrcieu.r-universe.dev | TwoSampleMR and other MRC-IEU packages | 200 | |
| files.pythonhosted.org, pypi.org | Python packages | 206 | Direct, no proxy |
| conda.anaconda.org | conda-forge / bioconda | 206 | |
| molgenis26.gcc.rug.nl | eQTLGen downloads | 200 | |
| storage.googleapis.com | FinnGen R10 public summary statistics | 200 | The manifest is `R10_manifest.tsv`. A first guess at the file name gave 404 from Google (host reachable) |
| api.semanticscholar.org | literature metadata | 200 | Unauthenticated, rate-limited |

## Consequences for the modules

- **Module 1 (MR) can run** without OpenGWAS:
  - outcome GWAS from the GWAS Catalog FTP or FinnGen R10;
  - eQTL exposures from eQTLGen;
  - gene coordinates from Ensembl;
  - TwoSampleMR from r-universe for the cross-check.
- **Literature:** open-access full text through the Europe PMC REST API and
  NCBI PMC. Publisher pages behind Cloudflare stay closed unless the
  researcher supplies the PDF.
