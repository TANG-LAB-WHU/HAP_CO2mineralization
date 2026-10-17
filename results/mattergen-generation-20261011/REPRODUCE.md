# Read-only verification and provenance

From this package directory, run the bundled standard-library verifier:

```console
python -B verify_results.py
```

It reads only this package. It checks manifest SHA-256 values, copied/transformed publication hashes, XYZ frame lengths and element counts against the recorded filters. It does not parse coordinates into geometry, calculate distances/energies/forces, import a model, query a GPU or network, or write output files. Passing it establishes package consistency, not scientific validity.

`source-provenance.json` uses paths relative to the original HAP_CO2mineralization project. Original byte hashes and transformations are explicit. The copied smoke manifests/structures, HAP44 filter/budget/structures and HAP44 screening contract are unmodified. Runtime/execution summaries whitelist relevant fields and omit GPU UUIDs, owner/container IDs, private absolute paths, host process inventories, mount commands, approval/import receipts and other operational details. Do not compare a transformed publication hash directly to its original-source hash.

Three original identifiers fix the later pilot:

- HAP44 filter summary: `69975a0b91a16dbc71cf7372707bb41c9a8eee944353f81ba28cd02e3300c6f1`.
- Screening contract: `8bb0668289c0c45ee30a5a6d5cb92e2627996d06073e090163d1a1b2d5c83369`.
- Full execution-parameter contract, retained locally but not copied here: `0d33bb1e6e7ce97b407cd8dd13295002900f200b66321f42e12e05da776ded5a`.

The running contracts under `.research-runs/campaign-002/mattergen-hap44/input/` matched the corresponding worktree files under `.research-runs/campaign-002/worktree/0.Generative_Discovery/`. The filters and raw candidate hashes are retained even though all candidates were rejected. The smoke filter's `structures_sha256` identifies its relaxed MatterSim structure; it is not the hash of the unrelaxed generated file. Other original manifest fields are preserved verbatim without reinterpreting their hash semantics.

Model checkpoints, third-party weights/reference datasets, original runtime source, large/native logs, telemetry and private approval records are excluded. Model and worker identities are hashes only; this bundle does not supply a runnable generation environment. Scientific reruns or new CSP training require a separate decision and their own reproducible configuration. Neither a frozen contract nor successful process termination proves scientific acceptance.
