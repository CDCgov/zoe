/label ~Release

- [ ] Make sure all needed PRs or issues are merged/closed for this cycle (consider what may be needed downstream)
- [ ] For each commit, ensure it has accurate documentation and a CHANGELOG entry if needed
- [ ] Run `cargo semver-checks --baseline-version <OLD VERSION> --release-type minor --all-features` to check for missing CHANGELOG entries
- [ ] Check for broken doc links in the pipeline output (and consider whether deny can be turned back on. See !463)
- [ ] Condense/reorder CHANGELOG
- [ ] Add date for the version being released in the CHANGELOG
- [ ] Make sure `Cargo.toml` version is correct
- [ ] Run `cargo update`
- [ ] Update `CITATION.bib`
- [ ] Get approval for release PR
- [ ] Squash/merge release PR
- [ ] Tag release
- [ ] Announce release in IFX Releases and Upgrades channel of ID Bioinformatics group on Teams
- [ ] Open PR updating to next `dev` version in TOML and CHANGELOG and removing any code deprecated before the release prior to this one