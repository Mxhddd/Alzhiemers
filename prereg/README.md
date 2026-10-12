# prereg/

Frozen pre-registrations. A hypothesis-testing module refuses to run on real
target data unless `python -m src.common.prereg check <name>` passes, which
requires:

1. `prereg/<name>.md` — the local copy of the pre-registration;
2. `prereg/<name>.freeze.yaml` — written by `make freeze-prereg NAME=<name> OSF=<url>`,
   holding the document's SHA-256, the UTC freeze time and the OSF URL;
3. the document unchanged since freezing (hash match);
4. an `https://osf.io/` URL — the third-party timestamp;
5. both files committed to git with no uncommitted changes.

To amend a pre-registration, write a new, separately named document (e.g.
`mediation_amendment1.md`), freeze it, and say in the write-up what changed and
when. Never edit a frozen document.
