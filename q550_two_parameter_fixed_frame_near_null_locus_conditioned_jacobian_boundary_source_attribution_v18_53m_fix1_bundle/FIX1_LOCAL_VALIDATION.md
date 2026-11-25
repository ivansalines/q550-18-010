# Q550 v18.53m-fix1 — Local validation, 2026-11-09

- Original v18.53m server archive intact: SHA-256 `a9348cf4eb3b9de75d75396dd3635e5ad68e2094272a5ad8b0e2e9c3b7b7c5dd`.
- Original v18.53m CRC passed, source manifest 16/16, output manifest 19/19; authentic v18.53l parent verified.
- Original failure is a *trace-domain exception*, not evidence for or against projected zeros: 1/10 Q550 IVPs and no trustworthy 2D root output.
- Fix1 source manifest verified. Original synthetic root/Jacobian/gauge test PASS; new trace-domain lower-bound and boundary-safe finite difference regression PASS.
- Full `bash scripts/run_on_server.sh --preflight` PASS; authentic v18.53l source manifest 15/15 and output manifest 28/28; no physical Q550 IVPs executed locally.
- Fixed origin range `[12.8,15.4]`, width range `[12.70,12.91]`; trace `[12.8,56.0]`, so endpoint maximum `54.13` remains inside coverage.
- This distribution ZIP deliberately excludes the generated local preflight report, which is not physical evidence. The server creates a fresh report and complete result ZIP at execution time.
