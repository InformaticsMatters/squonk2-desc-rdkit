# squonk2-desc-rdkit

Squonk2 Data Manager Jobs that calculate [RDKit] descriptors and fingerprints
for a set of molecules. Three Jobs share one image:

| Job | What it calculates |
| --- | ------------------ |
| `rdkit-descriptor-generator` | 2D RDKit descriptors |
| `rdkfp-generator` | RDKit fingerprints |
| `morgan-fingerprint-generator` | Morgan fingerprints |

The Job Definitions live in [`data-manager/jobs.yaml`](data-manager/jobs.yaml);
documentation for the Jobs is in [`data-manager/docs`](data-manager/docs).

## The image

The image is built from [`Dockerfile`](Dockerfile) and published as
`informaticsmatters/rdkit-descriptors`. Dependencies are managed with Poetry;
`poetry.lock` is what the image installs, so a dependency change means
relocking with the Poetry version the `Dockerfile` pins.

Note that RDKit uses CalVer, so a caret constraint does not mean what it
usually does: `^2025.9.1` expands to `>=2025.9.1,<2026.0.0` and silently
excludes every 2026 release. Prefer a plain `>=` floor.

## Testing

Jobs are tested with [jote]:

```bash
poetry install --only dev
jote
```

`jote --dry-run` validates the Job Definitions against the schema without
running anything, which is what CI does on every branch.

## Releasing

1. Run the `publish-tag` workflow with the new tag.
2. Set each Job Definition `version` and `image.tag` to that same value.
3. Cut the matching Git tag.

Never reuse a container tag — the Data Manager caches any tag other than
`latest`/`stable` per Kubernetes node. See `docs/versioning.md` in the
[squonk2-jobs] umbrella repository.

[jote]: https://github.com/InformaticsMatters/squonk2-data-manager-job-tester
[rdkit]: https://www.rdkit.org/
[squonk2-jobs]: https://github.com/InformaticsMatters/squonk2-jobs
