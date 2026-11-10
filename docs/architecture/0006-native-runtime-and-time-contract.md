# 0006 Native runtime capabilities and state-column timing

## Status

Accepted and implemented in the inspected development snapshot `3748e416d3a8120cc860119bd7224c6f9ad41ff7`. This record reconciles historical ADRs 0002–0004 with the native path; it does not retroactively declare all their modes-0–2 refactors complete.

## Context

Whole-member hooks remain useful for existing models. Large finite-element members need owned field contracts and persistent native contexts. File-backed analysis alone does not make the forecast distributed. Separately, initialization and observation/process-noise scheduling previously used inconsistent clocks. Recovery requires both state and policy identity.

## Decision

Retain legacy hooks and add native capabilities through the mode-3 registry. Use segmented field ownership, spatial model groups and matching ensemble communicators. Stream observation products and state rows through the existing selected transform. Preserve stable member identities. Expose inversion, restore and local analysis only when the application implements them.

Save the initial physical state at column 0, without forecasting. Forecast cycle k produces column k+1, applies scheduled process noise, and analyzes an observation attached to that column. Reject initial-time observations until a separate pre-forecast analysis stage exists. Persist `state-columns-v2` and a time-grid signature before state creation; refuse incompatible old histories rather than relabeling them.

Use sharded checkpoints, explicit manifests and bounded history publication. Preserve VDS backing files. On rank failure, report the original error and terminate without assuming unmatched cleanup collectives can succeed. Perform stopped-run recovery separately from failure handling.

## Consequences and limits

Mode 3 is directly selectable for a registered adapter. PIG has specialized native forecast/restart/history support; generic native applications have narrower restrictions. Infrastructure tests and one adapter's success do not certify all adapters or controlled scaling. The PIG noise bridge has a guarded global-size workspace, and the small Lorenz output can gather a tiny ensemble. Scientific error/factor, localization and inflation policies remain explicit experiment choices.

Use the current handbook and evidence registry for exact source, tests and measured limits. Release version/license metadata still require reconciliation before assigning a new public release. No new release number is introduced by this record.
