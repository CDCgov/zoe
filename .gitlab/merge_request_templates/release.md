/label ~release

- [ ] Make sure all needed PRs or issues are merged/closed for this cycle (consider what may be needed downstream)
- [ ] Make sure recent documentation is correct and complete
- [ ] Check for broken doc links in the pipeline output (and consider whether deny can be turned back on. See !463)
- [ ] Make sure latest CHANGELOG is correct
- [ ] Condense/reorder CHANGELOG
- [ ] Add date for the version being released in the CHANGELOG
- [ ] Make sure `Cargo.toml` version is correct
- [ ] Run `cargo update`
- [ ] Get approval for release PR
- [ ] Squash/merge release PR
- [ ] Tag release
- [ ] Announce release to colleagues if an important change has been made